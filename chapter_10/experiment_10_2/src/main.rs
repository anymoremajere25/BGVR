use std::collections::{HashMap, BTreeMap};
use std::path::PathBuf;
use std::fs::File;
use std::io::{Write, BufReader, BufRead};
use std::sync::Arc;
use anyhow::{Context, Result, bail};
use clap::Parser;
use log::{info, warn, error, debug};
use rayon::prelude::*;
use rust_htslib::bam::{self, Read, Record, HeaderView};
use serde::{Deserialize, Serialize};
use statrs::distribution::{Poisson, Continuous};
use itertools::Itertools;
use dashmap::DashMap;
use indicatif::{ProgressBar, ProgressStyle};
use chrono::{DateTime, Utc};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input BAM file path
    #[arg(short, long)]
    input: PathBuf,
    
    /// Output peaks file path (JSON format)
    #[arg(short, long, default_value = "peaks.json")]
    output: PathBuf,
    
    /// Optional BED file with predefined intervals to analyze
    #[arg(long)]
    intervals: Option<PathBuf>,
    
    /// Window size for genome-wide scanning (bp)
    #[arg(short, long, default_value_t = 500)]
    window_size: u64,
    
    /// Minimum coverage threshold for peak calling
    #[arg(short, long, default_value_t = 5.0)]
    min_coverage: f64,
    
    /// P-value threshold for peak significance
    #[arg(short, long, default_value_t = 0.05)]
    pvalue_threshold: f64,
    
    /// Fragment shift for ChIP-seq data (bp)
    #[arg(long, default_value_t = 75)]
    fragment_shift: i32,
    
    /// Extend reads to this length (bp)
    #[arg(long, default_value_t = 150)]
    extend_reads: u64,
    
    /// Number of threads to use
    #[arg(short, long, default_value_t = 4)]
    threads: usize,
    
    /// Enable verbose logging
    #[arg(short, long)]
    verbose: bool,
    
    /// Output BED format as well
    #[arg(long)]
    output_bed: bool,
    
    /// Cache intermediate results to disk
    #[arg(long)]
    cache_results: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
struct Interval {
    chrom: String,
    start: u64,
    end: u64,
    name: Option<String>,
}

impl Interval {
    fn new(chrom: String, start: u64, end: u64) -> Self {
        Self { chrom, start, end, name: None }
    }
    
    fn with_name(chrom: String, start: u64, end: u64, name: String) -> Self {
        Self { chrom, start, end, name: Some(name) }
    }
    
    fn overlaps(&self, other: &Interval) -> bool {
        self.chrom == other.chrom && self.start < other.end && self.end > other.start
    }
    
    fn length(&self) -> u64 {
        self.end.saturating_sub(self.start)
    }
}

#[derive(Debug, Serialize, Deserialize)]
struct Peak {
    interval: Interval,
    coverage: f64,
    normalized_coverage: f64,
    pvalue: f64,
    qvalue: f64,
    fold_enrichment: f64,
    summit: u64,
    summit_coverage: f64,
    read_count: u64,
    effective_length: u64,
}

#[derive(Debug)]
struct CoverageProfile {
    positions: BTreeMap<u64, f64>,
    total_reads: u64,
    mean_coverage: f64,
}

#[derive(Debug)]
struct IntervalTree {
    intervals_by_chrom: HashMap<String, Vec<Interval>>,
    index_by_chrom: HashMap<String, Vec<(u64, u64, usize)>>, // (start, end, index)
}

impl IntervalTree {
    fn new() -> Self {
        Self {
            intervals_by_chrom: HashMap::new(),
            index_by_chrom: HashMap::new(),
        }
    }
    
    fn insert(&mut self, interval: Interval) {
        let chrom = interval.chrom.clone();
        let intervals = self.intervals_by_chrom.entry(chrom.clone()).or_insert_with(Vec::new);
        let index = intervals.len();
        intervals.push(interval.clone());
        
        let index_vec = self.index_by_chrom.entry(chrom).or_insert_with(Vec::new);
        index_vec.push((interval.start, interval.end, index));
    }
    
    fn build_index(&mut self) {
        for (_, index_vec) in self.index_by_chrom.iter_mut() {
            index_vec.sort_by_key(|&(start, _, _)| start);
        }
    }
    
    fn query_overlaps(&self, query: &Interval) -> Vec<&Interval> {
        let mut overlaps = Vec::new();
        
        if let Some(intervals) = self.intervals_by_chrom.get(&query.chrom) {
            if let Some(index_vec) = self.index_by_chrom.get(&query.chrom) {
                for &(start, end, idx) in index_vec {
                    if start >= query.end {
                        break; // No more overlaps possible
                    }
                    if end > query.start {
                        overlaps.push(&intervals[idx]);
                    }
                }
            }
        }
        
        overlaps
    }
}

struct PeakCaller {
    args: Arc<Args>,
    interval_tree: Option<IntervalTree>,
    coverage_cache: Arc<DashMap<String, CoverageProfile>>,
}

impl PeakCaller {
    fn new(args: Args) -> Result<Self> {
        let mut interval_tree = None;
        
        // Load predefined intervals if provided
        if let Some(intervals_file) = &args.intervals {
            info!("Loading intervals from {:?}", intervals_file);
            interval_tree = Some(Self::load_intervals(intervals_file)?);
        }
        
        Ok(Self {
            args: Arc::new(args),
            interval_tree,
            coverage_cache: Arc::new(DashMap::new()),
        })
    }
    
    fn load_intervals(intervals_file: &PathBuf) -> Result<IntervalTree> {
        let mut tree = IntervalTree::new();
        let file = File::open(intervals_file)
            .with_context(|| format!("Failed to open intervals file: {:?}", intervals_file))?;
        
        let reader = BufReader::new(file);
        let mut count = 0;
        
        for line in reader.lines() {
            let line = line?;
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 3 {
                warn!("Skipping malformed line: {}", line);
                continue;
            }
            
            let chrom = fields[0].to_string();
            let start: u64 = fields[1].parse()
                .with_context(|| format!("Invalid start coordinate: {}", fields[1]))?;
            let end: u64 = fields[2].parse()
                .with_context(|| format!("Invalid end coordinate: {}", fields[2]))?;
            
            let interval = if fields.len() > 3 && !fields[3].is_empty() {
                Interval::with_name(chrom, start, end, fields[3].to_string())
            } else {
                Interval::new(chrom, start, end)
            };
            
            tree.insert(interval);
            count += 1;
        }
        
        tree.build_index();
        info!("Loaded {} intervals", count);
        Ok(tree)
    }
    
    fn generate_genome_windows(&self, bam_header: &HeaderView) -> Result<IntervalTree> {
        let mut tree = IntervalTree::new();
        let window_size = self.args.window_size;
        
        info!("Generating genome-wide windows of size {} bp", window_size);
        
        for i in 0..bam_header.target_count() {
            let chrom_name = String::from_utf8_lossy(bam_header.tid2name(i)).to_string();
            let chrom_len = bam_header.target_len(i).unwrap_or(0) as u64;
            
            for start in (0..chrom_len).step_by(window_size as usize) {
                let end = (start + window_size).min(chrom_len);
                let interval = Interval::with_name(
                    chrom_name.clone(),
                    start,
                    end,
                    format!("window_{}_{}", start, end)
                );
                tree.insert(interval);
            }
        }
        
        tree.build_index();
        Ok(tree)
    }
    
    fn compute_coverage_profile(&self, reader: &mut bam::Reader, header: &HeaderView, 
                               interval: &Interval) -> Result<CoverageProfile> {
        let mut positions = BTreeMap::new();
        let mut read_count = 0u64;
        
        // Set up region for fetching
        let tid = header.tid(interval.chrom.as_bytes());
        if tid.is_none() {
            return Ok(CoverageProfile {
                positions,
                total_reads: 0,
                mean_coverage: 0.0,
            });
        }
        
        let tid = tid.unwrap();
        
        // Fetch reads in the interval region
        if let Err(e) = reader.fetch((tid, interval.start as i64, interval.end as i64)) {
            warn!("Failed to fetch region {}:{}-{}: {}", 
                  interval.chrom, interval.start, interval.end, e);
            return Ok(CoverageProfile {
                positions,
                total_reads: 0,
                mean_coverage: 0.0,
            });
        }
        
        for result in reader.records() {
            let record = result?;
            
            // Skip unmapped, secondary, and duplicate reads
            if record.is_unmapped() || record.is_secondary() || record.is_duplicate() {
                continue;
            }
            
            // Apply fragment shift and extension
            let read_start = record.pos() as u64;
            let shifted_start = if record.is_reverse() {
                read_start.saturating_sub(self.args.fragment_shift as u64)
            } else {
                read_start + self.args.fragment_shift as u64
            };
            
            let read_end = shifted_start + self.args.extend_reads;
            
            // Check if read overlaps with interval
            if read_end >= interval.start && shifted_start <= interval.end {
                // Add coverage for each position
                let overlap_start = shifted_start.max(interval.start);
                let overlap_end = read_end.min(interval.end);
                
                for pos in overlap_start..overlap_end {
                    *positions.entry(pos).or_insert(0.0) += 1.0;
                }
                read_count += 1;
            }
        }
        
        let total_coverage: f64 = positions.values().sum();
        let mean_coverage = if interval.length() > 0 {
            total_coverage / interval.length() as f64
        } else {
            0.0
        };
        
        Ok(CoverageProfile {
            positions,
            total_reads: read_count,
            mean_coverage,
        })
    }
    
    fn call_peaks_from_intervals(&self, intervals: &IntervalTree, 
                                background_rate: f64) -> Result<Vec<Peak>> {
        info!("Calling peaks from {} chromosomes", intervals.intervals_by_chrom.len());
        
        let pb = ProgressBar::new(
            intervals.intervals_by_chrom.values().map(|v| v.len()).sum::<usize>() as u64
        );
        pb.set_style(ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos:>7}/{len:7} {msg}")
            .unwrap());
        
        let peaks: Result<Vec<Peak>> = intervals.intervals_by_chrom
            .par_iter()
            .flat_map(|(chrom, intervals)| {
                intervals.par_iter().map(|interval| {
                    pb.inc(1);
                    self.analyze_interval(interval, background_rate)
                }).collect::<Vec<_>>()
            })
            .collect();
        
        pb.finish_with_message("Peak calling completed");
        
        let mut all_peaks: Vec<Peak> = peaks?.into_iter().flatten().collect();
        
        // Apply multiple testing correction
        self.apply_multiple_testing_correction(&mut all_peaks);
        
        // Filter by significance and coverage
        all_peaks.retain(|peak| {
            peak.coverage >= self.args.min_coverage && 
            peak.qvalue <= self.args.pvalue_threshold
        });
        
        // Sort by significance
        all_peaks.sort_by(|a, b| a.qvalue.partial_cmp(&b.qvalue).unwrap_or(std::cmp::Ordering::Equal));
        
        Ok(all_peaks)
    }
    
    fn analyze_interval(&self, interval: &Interval, background_rate: f64) -> Result<Option<Peak>> {
        // Try to get from cache first
        let cache_key = format!("{}:{}:{}", interval.chrom, interval.start, interval.end);
        
        if let Some(cached) = self.coverage_cache.get(&cache_key) {
            return self.create_peak_from_profile(interval, &cached, background_rate);
        }
        
        // Compute coverage profile
        let mut reader = bam::Reader::from_path(&self.args.input)
            .with_context(|| format!("Failed to open BAM file: {:?}", self.args.input))?;
        
        let header = reader.header().clone();
        let profile = self.compute_coverage_profile(&mut reader, &header, interval)?;
        
        // Cache the result if enabled
        if self.args.cache_results {
            self.coverage_cache.insert(cache_key, profile);
        }
        
        let cached_profile = self.coverage_cache.get(&cache_key).unwrap();
        self.create_peak_from_profile(interval, &cached_profile, background_rate)
    }
    
    fn create_peak_from_profile(&self, interval: &Interval, 
                               profile: &CoverageProfile, 
                               background_rate: f64) -> Result<Option<Peak>> {
        if profile.total_reads == 0 || profile.mean_coverage < self.args.min_coverage {
            return Ok(None);
        }
        
        // Find summit (position with highest coverage)
        let (summit, summit_coverage) = profile.positions
            .iter()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(&pos, &cov)| (pos, cov))
            .unwrap_or((interval.start + interval.length() / 2, profile.mean_coverage));
        
        // Calculate statistics
        let poisson = Poisson::new(background_rate.max(0.1))
            .context("Failed to create Poisson distribution")?;
        
        let pvalue = 1.0 - poisson.cdf(profile.mean_coverage - 1.0);
        let fold_enrichment = profile.mean_coverage / background_rate.max(0.1);
        
        // Normalize coverage by interval length
        let normalized_coverage = profile.total_reads as f64 / interval.length() as f64 * 1000.0; // RPKM-like
        
        Ok(Some(Peak {
            interval: interval.clone(),
            coverage: profile.mean_coverage,
            normalized_coverage,
            pvalue,
            qvalue: pvalue, // Will be corrected later
            fold_enrichment,
            summit,
            summit_coverage,
            read_count: profile.total_reads,
            effective_length: interval.length(),
        }))
    }
    
    fn calculate_background_rate(&self) -> Result<f64> {
        info!("Calculating genome-wide background rate...");
        
        let mut reader = bam::Reader::from_path(&self.args.input)?;
        let header = reader.header().clone();
        
        let mut total_reads = 0u64;
        let mut total_length = 0u64;
        
        // Sample random regions to estimate background
        let sample_regions = 1000;
        let sample_size = 10000u64; // 10kb regions
        
        for i in 0..header.target_count().min(sample_regions / 100) {
            let chrom_len = header.target_len(i).unwrap_or(0) as u64;
            if chrom_len < sample_size {
                continue;
            }
            
            let chrom_name = String::from_utf8_lossy(header.tid2name(i)).to_string();
            
            // Sample multiple regions per chromosome
            for _ in 0..10 {
                let start = fastrand::u64(0..chrom_len.saturating_sub(sample_size));
                let end = start + sample_size;
                
                let interval = Interval::new(chrom_name.clone(), start, end);
                let profile = self.compute_coverage_profile(&mut reader, &header, &interval)?;
                
                total_reads += profile.total_reads;
                total_length += sample_size;
            }
        }
        
        let background_rate = if total_length > 0 {
            total_reads as f64 / total_length as f64
        } else {
            1.0 // Default fallback
        };
        
        info!("Estimated background rate: {:.4} reads per bp", background_rate);
        Ok(background_rate)
    }
    
    fn apply_multiple_testing_correction(&self, peaks: &mut [Peak]) {
        let n_tests = peaks.len() as f64;
        
        // Sort by p-value
        peaks.sort_by(|a, b| a.pvalue.partial_cmp(&b.pvalue).unwrap_or(std::cmp::Ordering::Equal));
        
        // Apply Benjamini-Hochberg correction
        for (i, peak) in peaks.iter_mut().enumerate() {
            let rank = (i + 1) as f64;
            peak.qvalue = (peak.pvalue * n_tests / rank).min(1.0);
        }
        
        // Enforce monotonicity
        for i in (0..peaks.len().saturating_sub(1)).rev() {
            if peaks[i].qvalue > peaks[i + 1].qvalue {
                peaks[i].qvalue = peaks[i + 1].qvalue;
            }
        }
    }
    
    fn write_peaks(&self, peaks: &[Peak]) -> Result<()> {
        info!("Writing {} peaks to {:?}", peaks.len(), self.args.output);
        
        // Write JSON output
        let output = serde_json::to_string_pretty(peaks)
            .context("Failed to serialize peaks to JSON")?;
        
        let mut file = File::create(&self.args.output)
            .with_context(|| format!("Failed to create output file: {:?}", self.args.output))?;
        
        file.write_all(output.as_bytes())
            .with_context(|| format!("Failed to write to output file: {:?}", self.args.output))?;
        
        // Write BED output if requested
        if self.args.output_bed {
            let bed_path = self.args.output.with_extension("bed");
            self.write_bed_format(peaks, &bed_path)?;
        }
        
        Ok(())
    }
    
    fn write_bed_format(&self, peaks: &[Peak], bed_path: &PathBuf) -> Result<()> {
        info!("Writing BED format to {:?}", bed_path);
        
        let mut file = File::create(bed_path)
            .with_context(|| format!("Failed to create BED file: {:?}", bed_path))?;
        
        // Write header
        writeln!(file, "track name=\"Peaks\" description=\"Peak calls\" itemRgb=\"On\"")?;
        
        for (i, peak) in peaks.iter().enumerate() {
            let score = ((-10.0 * peak.qvalue.log10()).min(1000.0).max(0.0)) as u32;
            let summit_offset = peak.summit.saturating_sub(peak.interval.start);
            
            writeln!(
                file,
                "{}\t{}\t{}\tpeak_{}\t{}\t.\t{:.2}\t{:.2e}\t{:.2e}\t{}",
                peak.interval.chrom,
                peak.interval.start,
                peak.interval.end,
                i + 1,
                score,
                peak.coverage,
                peak.pvalue,
                peak.qvalue,
                summit_offset
            )?;
        }
        
        Ok(())
    }
    
    pub fn run(&mut self) -> Result<()> {
        let start_time = Utc::now();
        info!("Starting peak calling analysis at {}", start_time.format("%Y-%m-%d %H:%M:%S"));
        
        // Validate input file
        if !self.args.input.exists() {
            bail!("Input BAM file does not exist: {:?}", self.args.input);
        }
        
        // Open BAM file and get header
        let mut reader = bam::Reader::from_path(&self.args.input)?;
        let header = reader.header().clone();
        
        info!("BAM file contains {} chromosomes", header.target_count());
        
        // Get intervals to analyze
        let intervals = if let Some(ref tree) = self.interval_tree {
            tree.clone()
        } else {
            self.generate_genome_windows(&header)?
        };
        
        // Calculate background rate
        let background_rate = self.calculate_background_rate()?;
        
        // Call peaks
        let peaks = self.call_peaks_from_intervals(&intervals, background_rate)?;
        
        // Write results
        self.write_peaks(&peaks)?;
        
        let end_time = Utc::now();
        let duration = end_time.signed_duration_since(start_time);
        
        info!("Peak calling completed in {:.1} seconds", duration.num_milliseconds() as f64 / 1000.0);
        info!("Found {} significant peaks", peaks.len());
        
        if !peaks.is_empty() {
            let avg_coverage: f64 = peaks.iter().map(|p| p.coverage).sum::<f64>() / peaks.len() as f64;
            let best_qvalue = peaks.iter().map(|p| p.qvalue).fold(1.0, f64::min);
            let max_fold_enrichment = peaks.iter().map(|p| p.fold_enrichment).fold(0.0, f64::max);
            
            info!("Summary statistics:");
            info!("  Average coverage: {:.2}", avg_coverage);
            info!("  Best q-value: {:.2e}", best_qvalue);
            info!("  Max fold enrichment: {:.2}", max_fold_enrichment);
        }
        
        Ok(())
    }
}

fn main() -> Result<()> {
    let args = Args::parse();
    
    // Initialize logging
    let log_level = if args.verbose { "debug" } else { "info" };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(log_level)).init();
    
    info!("Interval-based peak caller starting...");
    info!("Parameters:");
    info!("  Input: {:?}", args.input);
    info!("  Output: {:?}", args.output);
    info!("  Intervals file: {:?}", args.intervals.as_ref().unwrap_or(&PathBuf::from("genome-wide")));
    info!("  Window size: {} bp", args.window_size);
    info!("  Min coverage: {}", args.min_coverage);
    info!("  P-value threshold: {}", args.pvalue_threshold);
    info!("  Fragment shift: {} bp", args.fragment_shift);
    info!("  Read extension: {} bp", args.extend_reads);
    info!("  Threads: {}", args.threads);
    
    // Set up thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .context("Failed to initialize thread pool")?;
    
    // Run peak calling
    let mut peak_caller = PeakCaller::new(args)?;
    peak_caller.run()?;
    
    Ok(())
}