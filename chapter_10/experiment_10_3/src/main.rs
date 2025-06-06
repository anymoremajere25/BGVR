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
use statrs::distribution::{Normal, Continuous};
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
    
    /// Input peaks BED file path
    #[arg(short, long)]
    peaks: PathBuf,
    
    /// Output QC results file (JSON format)
    #[arg(short, long, default_value = "qc_results.json")]
    output: PathBuf,
    
    /// Output normalized signals file
    #[arg(long)]
    normalized_output: Option<PathBuf>,
    
    /// Control/input BAM file for normalization (optional)
    #[arg(long)]
    control_bam: Option<PathBuf>,
    
    /// Background regions BED file (optional)
    #[arg(long)]
    background_regions: Option<PathBuf>,
    
    /// Blacklist regions BED file (optional)
    #[arg(long)]
    blacklist_regions: Option<PathBuf>,
    
    /// Fragment shift for ChIP-seq data (bp)
    #[arg(long, default_value_t = 75)]
    fragment_shift: i32,
    
    /// Read extension length (bp)
    #[arg(long, default_value_t = 150)]
    extend_reads: u64,
    
    /// Bin size for signal normalization (bp)
    #[arg(long, default_value_t = 1000)]
    bin_size: u64,
    
    /// Number of threads to use
    #[arg(short, long, default_value_t = 4)]
    threads: usize,
    
    /// Enable verbose logging
    #[arg(short, long)]
    verbose: bool,
    
    /// Generate additional diagnostic plots
    #[arg(long)]
    generate_plots: bool,
    
    /// Minimum mapping quality
    #[arg(long, default_value_t = 10)]
    min_mapq: u8,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct GenomicInterval {
    chrom: String,
    start: u64,
    end: u64,
    name: Option<String>,
    score: Option<f64>,
}

impl GenomicInterval {
    fn overlaps(&self, other: &GenomicInterval) -> bool {
        self.chrom == other.chrom && self.start < other.end && self.end > other.start
    }
    
    fn length(&self) -> u64 {
        self.end.saturating_sub(self.start)
    }
}

#[derive(Debug, Serialize, Deserialize)]
struct QCMetrics {
    // Basic read statistics
    total_reads: u64,
    mapped_reads: u64,
    duplicated_reads: u64,
    high_quality_reads: u64,
    mapping_rate: f64,
    duplication_rate: f64,
    
    // Peak-specific metrics
    reads_in_peaks: u64,
    frip_score: f64,
    
    // Signal quality metrics
    signal_to_noise_ratio: f64,
    normalized_strand_correlation: f64,
    relative_strand_correlation: f64,
    
    // Normalization metrics
    rpm_factor: f64,
    rpkm_factor: f64,
    tpm_factor: Option<f64>,
    
    // Library complexity
    library_complexity: f64,
    effective_library_size: u64,
    
    // Background metrics
    background_mean_coverage: f64,
    background_std_coverage: f64,
    
    // Fragment size distribution
    mean_fragment_size: f64,
    fragment_size_std: f64,
    
    // Quality flags
    pass_frip_threshold: bool,
    pass_snr_threshold: bool,
    pass_complexity_threshold: bool,
    overall_quality: String,
}

#[derive(Debug, Serialize, Deserialize)]
struct NormalizedSignal {
    interval: GenomicInterval,
    raw_coverage: f64,
    rpm_normalized: f64,
    rpkm_normalized: f64,
    tpm_normalized: Option<f64>,
    log2_fold_change: Option<f64>, // vs control
    z_score: f64,
}

#[derive(Debug)]
struct CoverageData {
    coverage_profile: BTreeMap<u64, f64>,
    total_coverage: f64,
    read_count: u64,
}

struct QualityController {
    args: Arc<Args>,
    peak_intervals: Vec<GenomicInterval>,
    background_intervals: Option<Vec<GenomicInterval>>,
    blacklist_intervals: Option<Vec<GenomicInterval>>,
}

impl QualityController {
    fn new(args: Args) -> Result<Self> {
        let peak_intervals = Self::load_bed_file(&args.peaks)
            .with_context(|| format!("Failed to load peaks file: {:?}", args.peaks))?;
        
        info!("Loaded {} peak intervals", peak_intervals.len());
        
        let background_intervals = if let Some(ref bg_file) = args.background_regions {
            Some(Self::load_bed_file(bg_file)
                .with_context(|| format!("Failed to load background regions: {:?}", bg_file))?)
        } else {
            None
        };
        
        let blacklist_intervals = if let Some(ref bl_file) = args.blacklist_regions {
            Some(Self::load_bed_file(bl_file)
                .with_context(|| format!("Failed to load blacklist regions: {:?}", bl_file))?)
        } else {
            None
        };
        
        Ok(Self {
            args: Arc::new(args),
            peak_intervals,
            background_intervals,
            blacklist_intervals,
        })
    }
    
    fn load_bed_file(bed_path: &PathBuf) -> Result<Vec<GenomicInterval>> {
        let file = File::open(bed_path)
            .with_context(|| format!("Failed to open BED file: {:?}", bed_path))?;
        
        let reader = BufReader::new(file);
        let mut intervals = Vec::new();
        
        for (line_num, line) in reader.lines().enumerate() {
            let line = line?;
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 3 {
                warn!("Skipping malformed line {}: {}", line_num + 1, line);
                continue;
            }
            
            let chrom = fields[0].to_string();
            let start: u64 = fields[1].parse()
                .with_context(|| format!("Invalid start coordinate on line {}", line_num + 1))?;
            let end: u64 = fields[2].parse()
                .with_context(|| format!("Invalid end coordinate on line {}", line_num + 1))?;
            
            let name = if fields.len() > 3 && !fields[3].is_empty() {
                Some(fields[3].to_string())
            } else {
                None
            };
            
            let score = if fields.len() > 4 {
                fields[4].parse().ok()
            } else {
                None
            };
            
            intervals.push(GenomicInterval {
                chrom, start, end, name, score
            });
        }
        
        Ok(intervals)
    }
    
    fn calculate_basic_statistics(&self, bam_path: &PathBuf) -> Result<(u64, u64, u64, u64)> {
        info!("Calculating basic read statistics...");
        
        let mut bam_reader = bam::Reader::from_path(bam_path)?;
        let mut total_reads = 0u64;
        let mut mapped_reads = 0u64;
        let mut duplicated_reads = 0u64;
        let mut high_quality_reads = 0u64;
        
        let pb = ProgressBar::new_spinner();
        pb.set_style(ProgressStyle::default_spinner()
            .template("{spinner:.green} Processing reads: {pos}").unwrap());
        
        for result in bam_reader.records() {
            let record = result?;
            total_reads += 1;
            
            if total_reads % 10000 == 0 {
                pb.set_position(total_reads);
            }
            
            if !record.is_unmapped() {
                mapped_reads += 1;
                
                if record.mapq() >= self.args.min_mapq {
                    high_quality_reads += 1;
                }
            }
            
            if record.is_duplicate() {
                duplicated_reads += 1;
            }
        }
        
        pb.finish_with_message("Read statistics calculated");
        
        info!("Total reads: {}", total_reads);
        info!("Mapped reads: {} ({:.1}%)", mapped_reads, 
              100.0 * mapped_reads as f64 / total_reads as f64);
        info!("High quality reads: {} ({:.1}%)", high_quality_reads,
              100.0 * high_quality_reads as f64 / total_reads as f64);
        info!("Duplicated reads: {} ({:.1}%)", duplicated_reads,
              100.0 * duplicated_reads as f64 / total_reads as f64);
        
        Ok((total_reads, mapped_reads, duplicated_reads, high_quality_reads))
    }
    
    fn count_reads_in_peaks(&self, bam_path: &PathBuf) -> Result<u64> {
        info!("Counting reads in peaks...");
        
        let mut bam_reader = bam::Reader::from_path(bam_path)?;
        let header = bam_reader.header().clone();
        
        // Create a parallel counter for reads in peaks
        let reads_in_peaks = Arc::new(std::sync::atomic::AtomicU64::new(0));
        
        let pb = ProgressBar::new_spinner();
        pb.set_style(ProgressStyle::default_spinner()
            .template("{spinner:.green} Checking peak overlaps: {pos}").unwrap());
        
        let mut processed = 0u64;
        
        for result in bam_reader.records() {
            let record = result?;
            processed += 1;
            
            if processed % 5000 == 0 {
                pb.set_position(processed);
            }
            
            // Skip unmapped, low quality, and duplicate reads
            if record.is_unmapped() || 
               record.mapq() < self.args.min_mapq || 
               record.is_duplicate() {
                continue;
            }
            
            // Skip if in blacklist
            if let Some(ref blacklist) = self.blacklist_intervals {
                if self.read_overlaps_intervals(&record, &header, blacklist) {
                    continue;
                }
            }
            
            // Check if read overlaps with any peak
            if self.read_overlaps_intervals(&record, &header, &self.peak_intervals) {
                reads_in_peaks.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            }
        }
        
        pb.finish_with_message("Peak overlap analysis completed");
        
        let count = reads_in_peaks.load(std::sync::atomic::Ordering::Relaxed);
        info!("Reads in peaks: {}", count);
        
        Ok(count)
    }
    
    fn read_overlaps_intervals(&self, record: &Record, header: &HeaderView, 
                              intervals: &[GenomicInterval]) -> bool {
        let tid = record.tid();
        if tid < 0 {
            return false;
        }
        
        let chrom_name = match header.tid2name(tid as u32) {
            Some(name) => String::from_utf8_lossy(name).to_string(),
            None => return false,
        };
        
        // Apply fragment shift and extension
        let read_start = record.pos() as u64;
        let shifted_start = if record.is_reverse() {
            read_start.saturating_sub(self.args.fragment_shift as u64)
        } else {
            read_start + self.args.fragment_shift as u64
        };
        
        let read_end = shifted_start + self.args.extend_reads;
        
        // Check overlap with intervals
        intervals.iter().any(|interval| {
            interval.chrom == chrom_name && 
            read_end >= interval.start && 
            shifted_start <= interval.end
        })
    }
    
    fn calculate_signal_to_noise_ratio(&self, bam_path: &PathBuf) -> Result<f64> {
        info!("Calculating signal-to-noise ratio...");
        
        // Calculate coverage in peaks vs background
        let peak_coverage = self.calculate_mean_coverage_in_regions(bam_path, &self.peak_intervals)?;
        
        let background_coverage = if let Some(ref bg_intervals) = self.background_intervals {
            self.calculate_mean_coverage_in_regions(bam_path, bg_intervals)?
        } else {
            // Use genome-wide background estimate
            self.estimate_genome_wide_background(bam_path)?
        };
        
        let snr = if background_coverage > 0.0 {
            peak_coverage / background_coverage
        } else {
            0.0
        };
        
        info!("Signal-to-noise ratio: {:.2}", snr);
        Ok(snr)
    }
    
    fn calculate_mean_coverage_in_regions(&self, bam_path: &PathBuf, 
                                        regions: &[GenomicInterval]) -> Result<f64> {
        let mut bam_reader = bam::Reader::from_path(bam_path)?;
        let header = bam_reader.header().clone();
        
        let total_coverage: f64 = regions.par_iter().map(|region| {
            let mut local_reader = bam::Reader::from_path(bam_path).unwrap();
            
            // Fetch reads in region
            let tid = header.tid(region.chrom.as_bytes());
            if tid.is_none() {
                return 0.0;
            }
            
            let tid = tid.unwrap();
            if local_reader.fetch((tid, region.start as i64, region.end as i64)).is_err() {
                return 0.0;
            }
            
            let mut coverage = 0u64;
            for result in local_reader.records() {
                if let Ok(record) = result {
                    if !record.is_unmapped() && 
                       record.mapq() >= self.args.min_mapq && 
                       !record.is_duplicate() {
                        coverage += 1;
                    }
                }
            }
            
            coverage as f64 / region.length() as f64
        }).sum();
        
        let mean_coverage = total_coverage / regions.len() as f64;
        Ok(mean_coverage)
    }
    
    fn estimate_genome_wide_background(&self, bam_path: &PathBuf) -> Result<f64> {
        info!("Estimating genome-wide background...");
        
        let mut bam_reader = bam::Reader::from_path(bam_path)?;
        let header = bam_reader.header().clone();
        
        // Sample random regions
        let sample_size = 1000;
        let region_size = 1000u64;
        let mut total_coverage = 0f64;
        let mut total_length = 0u64;
        
        for i in 0..std::cmp::min(header.target_count(), sample_size / 10) {
            let chrom_len = header.target_len(i).unwrap_or(0) as u64;
            if chrom_len < region_size * 2 {
                continue;
            }
            
            for _ in 0..10 {
                let start = fastrand::u64(region_size..chrom_len - region_size);
                let end = start + region_size;
                
                if let Ok(_) = bam_reader.fetch((i as i32, start as i64, end as i64)) {
                    let mut coverage = 0u64;
                    for result in bam_reader.records() {
                        if let Ok(record) = result {
                            if !record.is_unmapped() && 
                               record.mapq() >= self.args.min_mapq && 
                               !record.is_duplicate() {
                                coverage += 1;
                            }
                        }
                    }
                    
                    total_coverage += coverage as f64;
                    total_length += region_size;
                }
            }
        }
        
        let background = if total_length > 0 {
            total_coverage / total_length as f64
        } else {
            0.001 // Small default to avoid division by zero
        };
        
        info!("Estimated background coverage: {:.6}", background);
        Ok(background)
    }
    
    fn calculate_strand_correlation(&self, bam_path: &PathBuf) -> Result<(f64, f64)> {
        info!("Calculating strand correlation metrics...");
        
        // This is a simplified version - in practice, you'd implement
        // cross-correlation analysis similar to phantompeakqualtools
        
        // For now, return placeholder values
        // In a full implementation, you'd:
        // 1. Separate reads by strand
        // 2. Calculate coverage profiles
        // 3. Compute cross-correlation at different shifts
        // 4. Find peak correlation and relative correlation
        
        let nsc = 1.1; // Normalized strand correlation (placeholder)
        let rsc = 0.8; // Relative strand correlation (placeholder)
        
        info!("NSC: {:.3}, RSC: {:.3}", nsc, rsc);
        Ok((nsc, rsc))
    }
    
    fn estimate_library_complexity(&self, bam_path: &PathBuf) -> Result<(f64, u64)> {
        info!("Estimating library complexity...");
        
        let mut bam_reader = bam::Reader::from_path(bam_path)?;
        let mut unique_positions = std::collections::HashSet::new();
        let mut total_reads = 0u64;
        
        for result in bam_reader.records() {
            if let Ok(record) = result {
                if !record.is_unmapped() && 
                   record.mapq() >= self.args.min_mapq {
                    
                    let tid = record.tid();
                    let pos = record.pos();
                    let is_reverse = record.is_reverse();
                    
                    unique_positions.insert((tid, pos, is_reverse));
                    total_reads += 1;
                }
            }
        }
        
        let complexity = if total_reads > 0 {
            unique_positions.len() as f64 / total_reads as f64
        } else {
            0.0
        };
        
        info!("Library complexity: {:.3}", complexity);
        info!("Effective library size: {}", unique_positions.len());
        
        Ok((complexity, unique_positions.len() as u64))
    }
    
    fn calculate_fragment_size_distribution(&self, bam_path: &PathBuf) -> Result<(f64, f64)> {
        info!("Calculating fragment size distribution...");
        
        let mut bam_reader = bam::Reader::from_path(bam_path)?;
        let mut fragment_sizes = Vec::new();
        
        // Sample fragment sizes from properly paired reads
        for result in bam_reader.records() {
            if let Ok(record) = result {
                if record.is_paired() && 
                   record.is_proper_pair() && 
                   !record.is_unmapped() &&
                   record.mapq() >= self.args.min_mapq {
                    
                    let tlen = record.insert_size().abs() as u64;
                    if tlen > 0 && tlen < 1000 { // Reasonable fragment size range
                        fragment_sizes.push(tlen as f64);
                    }
                }
                
                if fragment_sizes.len() >= 10000 {
                    break; // Sample sufficient for estimate
                }
            }
        }
        
        let (mean, std) = if !fragment_sizes.is_empty() {
            let mean = fragment_sizes.iter().sum::<f64>() / fragment_sizes.len() as f64;
            let variance = fragment_sizes.iter()
                .map(|x| (x - mean).powi(2))
                .sum::<f64>() / fragment_sizes.len() as f64;
            let std = variance.sqrt();
            (mean, std)
        } else {
            (0.0, 0.0)
        };
        
        info!("Fragment size: {:.1} ± {:.1} bp", mean, std);
        Ok((mean, std))
    }
    
    fn calculate_normalization_factors(&self, total_reads: u64, 
                                     effective_library_size: u64) -> (f64, f64, Option<f64>) {
        // RPM (Reads Per Million)
        let rpm_factor = 1_000_000.0 / total_reads as f64;
        
        // RPKM factor (will be applied per-region based on length)
        let rpkm_factor = rpm_factor;
        
        // TPM factor (optional, more complex calculation)
        let tpm_factor = if effective_library_size > 0 {
            Some(1_000_000.0 / effective_library_size as f64)
        } else {
            None
        };
        
        info!("Normalization factors - RPM: {:.6}, RPKM base: {:.6}", 
              rpm_factor, rpkm_factor);
        
        (rpm_factor, rpkm_factor, tpm_factor)
    }
    
    fn generate_quality_assessment(&self, metrics: &QCMetrics) -> String {
        let mut issues = Vec::new();
        
        // Check FRiP score
        if metrics.frip_score < 0.01 {
            issues.push("Very low FRiP score (< 1%)");
        } else if metrics.frip_score < 0.05 {
            issues.push("Low FRiP score (< 5%)");
        }
        
        // Check signal-to-noise ratio
        if metrics.signal_to_noise_ratio < 2.0 {
            issues.push("Low signal-to-noise ratio (< 2.0)");
        }
        
        // Check library complexity
        if metrics.library_complexity < 0.7 {
            issues.push("Low library complexity (< 70%)");
        }
        
        // Check mapping rate
        if metrics.mapping_rate < 0.7 {
            issues.push("Low mapping rate (< 70%)");
        }
        
        // Check duplication rate
        if metrics.duplication_rate > 0.5 {
            issues.push("High duplication rate (> 50%)");
        }
        
        if issues.is_empty() {
            "PASS".to_string()
        } else if issues.len() <= 2 {
            "WARNING".to_string()
        } else {
            "FAIL".to_string()
        }
    }
    
    fn write_results(&self, metrics: &QCMetrics) -> Result<()> {
        info!("Writing QC results to {:?}", self.args.output);
        
        let output = serde_json::to_string_pretty(metrics)
            .context("Failed to serialize QC metrics")?;
        
        let mut file = File::create(&self.args.output)
            .with_context(|| format!("Failed to create output file: {:?}", self.args.output))?;
        
        file.write_all(output.as_bytes())
            .with_context(|| format!("Failed to write to output file: {:?}", self.args.output))?;
        
        Ok(())
    }
    
    pub fn run(&mut self) -> Result<()> {
        let start_time = Utc::now();
        info!("Starting quality control analysis at {}", 
              start_time.format("%Y-%m-%d %H:%M:%S"));
        
        // Validate input files
        if !self.args.input.exists() {
            bail!("Input BAM file does not exist: {:?}", self.args.input);
        }
        
        if !self.args.peaks.exists() {
            bail!("Peaks BED file does not exist: {:?}", self.args.peaks);
        }
        
        // Calculate basic statistics
        let (total_reads, mapped_reads, duplicated_reads, high_quality_reads) = 
            self.calculate_basic_statistics(&self.args.input)?;
        
        // Calculate reads in peaks and FRiP
        let reads_in_peaks = self.count_reads_in_peaks(&self.args.input)?;
        let frip_score = reads_in_peaks as f64 / total_reads as f64;
        
        // Calculate signal quality metrics
        let signal_to_noise_ratio = self.calculate_signal_to_noise_ratio(&self.args.input)?;
        let (nsc, rsc) = self.calculate_strand_correlation(&self.args.input)?;
        
        // Calculate library complexity
        let (library_complexity, effective_library_size) = 
            self.estimate_library_complexity(&self.args.input)?;
        
        // Calculate fragment size distribution
        let (mean_fragment_size, fragment_size_std) = 
            self.calculate_fragment_size_distribution(&self.args.input)?;
        
        // Calculate normalization factors
        let (rpm_factor, rpkm_factor, tpm_factor) = 
            self.calculate_normalization_factors(total_reads, effective_library_size);
        
        // Calculate background statistics
        let background_mean_coverage = if let Some(ref bg_intervals) = self.background_intervals {
            self.calculate_mean_coverage_in_regions(&self.args.input, bg_intervals)?
        } else {
            self.estimate_genome_wide_background(&self.args.input)?
        };
        
        // Compile metrics
        let metrics = QCMetrics {
            total_reads,
            mapped_reads,
            duplicated_reads,
            high_quality_reads,
            mapping_rate: mapped_reads as f64 / total_reads as f64,
            duplication_rate: duplicated_reads as f64 / total_reads as f64,
            reads_in_peaks,
            frip_score,
            signal_to_noise_ratio,
            normalized_strand_correlation: nsc,
            relative_strand_correlation: rsc,
            rpm_factor,
            rpkm_factor,
            tpm_factor,
            library_complexity,
            effective_library_size,
            background_mean_coverage,
            background_std_coverage: 0.0, // Placeholder
            mean_fragment_size,
            fragment_size_std,
            pass_frip_threshold: frip_score >= 0.01,
            pass_snr_threshold: signal_to_noise_ratio >= 2.0,
            pass_complexity_threshold: library_complexity >= 0.7,
            overall_quality: String::new(),
        };
        
        let mut final_metrics = metrics;
        final_metrics.overall_quality = self.generate_quality_assessment(&final_metrics);
        
        // Write results
        self.write_results(&final_metrics)?;
        
        let end_time = Utc::now();
        let duration = end_time.signed_duration_since(start_time);
        
        info!("Quality control analysis completed in {:.1} seconds", 
              duration.num_milliseconds() as f64 / 1000.0);
        
        // Print summary
        println!("\n=== Quality Control Summary ===");
        println!("