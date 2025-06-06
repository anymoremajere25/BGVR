use std::collections::HashMap;
use std::path::PathBuf;
use std::fs::File;
use std::io::{Write, BufReader, BufRead};
use anyhow::{Context, Result, bail};
use clap::Parser;
use log::{info, warn};
use rayon::prelude::*;
use rust_htslib::bam::{self, Read, Record, HeaderView};
use serde::{Deserialize, Serialize};
use statrs::distribution::{Poisson, Continuous};

#[derive(Parser, Debug)]
#[command(author, version, about = "Shift-aware peak caller for ChIP-seq data")]
struct Args {
    /// Input BAM file
    bam_file: PathBuf,
    
    /// Shift estimate JSON file
    shift_file: PathBuf,
    
    /// Output peaks BED file
    #[arg(short, long, default_value = "peaks.bed")]
    output: PathBuf,
    
    /// Window size for peak calling (bp)
    #[arg(short, long, default_value_t = 200)]
    window_size: u64,
    
    /// Minimum coverage threshold
    #[arg(short, long, default_value_t = 5.0)]
    min_coverage: f64,
    
    /// P-value threshold for significance
    #[arg(short, long, default_value_t = 0.05)]
    pvalue_threshold: f64,
    
    /// Minimum peak width (bp)
    #[arg(long, default_value_t = 50)]
    min_peak_width: u64,
    
    /// Maximum peak width (bp)
    #[arg(long, default_value_t = 5000)]
    max_peak_width: u64,
    
    /// Minimum mapping quality
    #[arg(long, default_value_t = 10)]
    min_mapq: u8,
    
    /// Number of threads
    #[arg(short, long, default_value_t = 4)]
    threads: usize,
    
    /// Enable verbose logging
    #[arg(short, long)]
    verbose: bool,
}

#[derive(Debug, Deserialize)]
struct ShiftEstimate {
    estimated_shift: i64,
    confidence_score: f64,
}

#[derive(Debug, Serialize)]
struct Peak {
    chrom: String,
    start: u64,
    end: u64,
    name: String,
    score: u32,
    strand: char,
    signal_value: f64,
    pvalue: f64,
    qvalue: f64,
    summit_offset: u64,
}

struct ShiftAwarePeakCaller {
    args: Args,
    fragment_shift: i64,
    confidence_score: f64,
}

impl ShiftAwarePeakCaller {
    fn new(args: Args) -> Result<Self> {
        // Load shift estimate
        let shift_data = Self::load_shift_estimate(&args.shift_file)?;
        
        info!("Loaded fragment shift: {} bp (confidence: {:.3})", 
              shift_data.estimated_shift, shift_data.confidence_score);
        
        if shift_data.confidence_score < 0.5 {
            warn!("Low confidence shift estimate ({}). Results may be unreliable.", 
                  shift_data.confidence_score);
        }
        
        Ok(Self {
            args,
            fragment_shift: shift_data.estimated_shift,
            confidence_score: shift_data.confidence_score,
        })
    }
    
    fn load_shift_estimate(shift_file: &PathBuf) -> Result<ShiftEstimate> {
        let file = File::open(shift_file)
            .with_context(|| format!("Failed to open shift file: {:?}", shift_file))?;
        
        let reader = BufReader::new(file);
        let shift_data: ShiftEstimate = serde_json::from_reader(reader)
            .with_context(|| format!("Failed to parse shift file: {:?}", shift_file))?;
        
        Ok(shift_data)
    }
    
    fn calculate_shifted_coverage(&self) -> Result<HashMap<String, Vec<(u64, f64)>>> {
        info!("Calculating coverage with fragment shift correction...");
        
        let mut bam_reader = bam::Reader::from_path(&self.args.bam_file)?;
        let header = bam_reader.header().clone();
        
        // Coverage maps for each chromosome
        let mut coverage_maps: HashMap<String, HashMap<u64, u64>> = HashMap::new();
        
        let mut reads_processed = 0u64;
        let mut reads_shifted = 0u64;
        
        for result in bam_reader.records() {
            let record = result?;
            reads_processed += 1;
            
            if reads_processed % 100000 == 0 {
                info!("Processed {} reads, {} shifted", reads_processed, reads_shifted);
            }
            
            // Filter reads
            if record.is_unmapped() || 
               record.is_secondary() || 
               record.is_duplicate() ||
               record.mapq() < self.args.min_mapq {
                continue;
            }
            
            // Get chromosome name
            let tid = record.tid();
            if tid < 0 {
                continue;
            }
            
            let chrom_name = match header.tid2name(tid as u32) {
                Some(name) => String::from_utf8_lossy(name).to_string(),
                None => continue,
            };
            
            // Apply fragment shift
            let original_pos = record.pos() as u64;
            let shifted_pos = if record.is_reverse() {
                // For reverse strand, shift in opposite direction
                original_pos.saturating_sub(self.fragment_shift as u64)
            } else {
                // For forward strand, shift towards 3' end
                original_pos + self.fragment_shift as u64
            };
            
            reads_shifted += 1;
            
            // Calculate window
            let window_start = (shifted_pos / self.args.window_size) * self.args.window_size;
            
            // Update coverage
            let chrom_coverage = coverage_maps.entry(chrom_name).or_insert_with(HashMap::new);
            *chrom_coverage.entry(window_start).or_insert(0) += 1;
        }
        
        info!("Processed {} reads total, {} passed filters and were shifted", 
              reads_processed, reads_shifted);
        
        // Convert to sorted vectors
        let coverage_data: HashMap<String, Vec<(u64, f64)>> = coverage_maps
            .into_iter()
            .map(|(chrom, windows)| {
                let mut sorted_windows: Vec<(u64, f64)> = windows
                    .into_iter()
                    .map(|(pos, count)| (pos, count as f64))
                    .collect();
                sorted_windows.sort_by_key(|&(pos, _)| pos);
                (chrom, sorted_windows)
            })
            .collect();
        
        Ok(coverage_data)
    }
    
    fn call_peaks(&self, coverage_data: HashMap<String, Vec<(u64, f64)>>) -> Result<Vec<Peak>> {
        info!("Calling peaks with statistical significance testing...");
        
        let mut all_peaks = Vec::new();
        
        for (chrom, windows) in coverage_data {
            if windows.is_empty() {
                continue;
            }
            
            // Calculate background rate for this chromosome
            let total_coverage: f64 = windows.iter().map(|(_, cov)| cov).sum();
            let background_rate = total_coverage / windows.len() as f64;
            
            if background_rate < 0.1 {
                continue;
            }
            
            info!("Processing chromosome {} with background rate: {:.2}", 
                  chrom, background_rate);
            
            // Use Poisson distribution for significance testing
            let poisson = Poisson::new(background_rate)
                .context("Failed to create Poisson distribution")?;
            
            // Find significant windows
            let mut significant_windows = Vec::new();
            for (pos, coverage) in &windows {
                if *coverage >= self.args.min_coverage {
                    let pvalue = 1.0 - poisson.cdf(*coverage - 1.0);
                    if pvalue <= self.args.pvalue_threshold {
                        significant_windows.push((*pos, *coverage, pvalue));
                    }
                }
            }
            
            // Merge nearby significant windows into peaks
            let chromosome_peaks = self.merge_windows_into_peaks(
                chrom.clone(), 
                significant_windows, 
                background_rate
            );
            
            all_peaks.extend(chromosome_peaks);
        }
        
        // Sort peaks by significance
        all_peaks.sort_by(|a, b| a.pvalue.partial_cmp(&b.pvalue).unwrap_or(std::cmp::Ordering::Equal));
        
        // Apply multiple testing correction (Benjamini-Hochberg)
        let n_tests = all_peaks.len() as f64;
        for (i, peak) in all_peaks.iter_mut().enumerate() {
            let rank = (i + 1) as f64;
            peak.qvalue = (peak.pvalue * n_tests / rank).min(1.0);
        }
        
        // Enforce monotonicity in q-values
        for i in (0..all_peaks.len().saturating_sub(1)).rev() {
            if all_peaks[i].qvalue > all_peaks[i + 1].qvalue {
                all_peaks[i].qvalue = all_peaks[i + 1].qvalue;
            }
        }
        
        info!("Found {} significant peaks", all_peaks.len());
        Ok(all_peaks)
    }
    
    fn merge_windows_into_peaks(&self, chrom: String, 
                               significant_windows: Vec<(u64, f64, f64)>,
                               background_rate: f64) -> Vec<Peak> {
        if significant_windows.is_empty() {
            return Vec::new();
        }
        
        let mut peaks = Vec::new();
        let mut current_peak_windows = vec![significant_windows[0]];
        
        for &(pos, coverage, pvalue) in significant_windows.iter().skip(1) {
            let last_pos = current_peak_windows.last().unwrap().0;
            
            // Check if this window is adjacent or nearby
            if pos <= last_pos + 2 * self.args.window_size {
                // Extend current peak
                current_peak_windows.push((pos, coverage, pvalue));
            } else {
                // Finalize current peak and start new one
                if let Some(peak) = self.finalize_peak(
                    chrom.clone(), 
                    &current_peak_windows, 
                    background_rate, 
                    peaks.len() + 1
                ) {
                    peaks.push(peak);
                }
                current_peak_windows = vec![(pos, coverage, pvalue)];
            }
        }
        
        // Finalize last peak
        if let Some(peak) = self.finalize_peak(
            chrom.clone(), 
            &current_peak_windows, 
            background_rate, 
            peaks.len() + 1
        ) {
            peaks.push(peak);
        }
        
        peaks
    }
    
    fn finalize_peak(&self, chrom: String, windows: &[(u64, f64, f64)], 
                    background_rate: f64, peak_number: usize) -> Option<Peak> {
        if windows.is_empty() {
            return None;
        }
        
        let start = windows.first().unwrap().0;
        let end = windows.last().unwrap().0 + self.args.window_size;
        let peak_width = end - start;
        
        // Filter by peak width
        if peak_width < self.args.min_peak_width || peak_width > self.args.max_peak_width {
            return None;
        }
        
        // Find summit (window with highest coverage)
        let (summit_pos, max_coverage, best_pvalue) = windows
            .iter()
            .max_by(|(_, a, _), (_, b, _)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .map(|&(pos, cov, pval)| (pos, cov, pval))?;
        
        let summit_offset = summit_pos - start;
        
        // Calculate aggregated statistics
        let total_coverage: f64 = windows.iter().map(|(_, cov, _)| cov).sum();
        let mean_coverage = total_coverage / windows.len() as f64;
        
        // Use the most significant p-value
        let peak_pvalue = windows.iter()
            .map(|(_, _, pval)| pval)
            .fold(1.0, |a, &b| a.min(b));
        
        // Calculate score (capped at 1000 for BED format)
        let score = ((-10.0 * peak_pvalue.log10()).min(1000.0).max(0.0)) as u32;
        
        Some(Peak {
            chrom,
            start,
            end,
            name: format!("peak_{}", peak_number),
            score,
            strand: '.',
            signal_value: mean_coverage,
            pvalue: peak_pvalue,
            qvalue: peak_pvalue, // Will be corrected later
            summit_offset,
        })
    }
    
    fn write_peaks(&self, peaks: &[Peak]) -> Result<()> {
        info!("Writing {} peaks to {:?}", peaks.len(), self.args.output);
        
        let mut file = File::create(&self.args.output)
            .with_context(|| format!("Failed to create output file: {:?}", self.args.output))?;
        
        // Write header comment
        writeln!(file, "# Peaks called with fragment shift: {} bp (confidence: {:.3})", 
                self.fragment_shift, self.confidence_score)?;
        writeln!(file, "# chrom\tstart\tend\tname\tscore\tstrand\tsignal\tpValue\tqValue\tsummit")?;
        
        for peak in peaks {
            writeln!(
                file,
                "{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2e}\t{:.2e}\t{}",
                peak.chrom,
                peak.start,
                peak.end,
                peak.name,
                peak.score,
                peak.strand,
                peak.signal_value,
                peak.pvalue,
                peak.qvalue,
                peak.summit_offset
            )?;
        }
        
        Ok(())
    }
    
    pub fn run(&mut self) -> Result<()> {
        info!("Starting shift-aware peak calling...");
        
        // Validate input files
        if !self.args.bam_file.exists() {
            bail!("BAM file does not exist: {:?}", self.args.bam_file);
        }
        
        // Calculate coverage with shift correction
        let coverage_data = self.calculate_shifted_coverage()?;
        
        // Call peaks
        let peaks = self.call_peaks(coverage_data)?;
        
        // Write results
        self.write_peaks(&peaks)?;
        
        // Print summary
        println!("\n=== Peak Calling Summary ===");
        println!("Fragment shift applied: {} bp", self.fragment_shift);
        println!("Shift confidence: {:.3}", self.confidence_score);
        println!("Peaks called: {}", peaks.len());
        if !peaks.is_empty() {
            let avg_signal: f64 = peaks.iter().map(|p| p.signal_value).sum::<f64>() / peaks.len() as f64;
            let best_pvalue = peaks.iter().map(|p| p.pvalue).fold(1.0, f64::min);
            println!("Average signal: {:.2}", avg_signal);
            println!("Best p-value: {:.2e}", best_pvalue);
        }
        println!("===========================");
        
        Ok(())
    }
}

fn main() -> Result<()> {
    let args = Args::parse();
    
    // Initialize logging
    let log_level = if args.verbose { "debug" } else { "info" };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(log_level)).init();
    
    info!("Shift-aware peak caller starting...");
    
    // Set thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .context("Failed to initialize thread pool")?;
    
    // Run peak calling
    let mut peak_caller = ShiftAwarePeakCaller::new(args)?;
    peak_caller.run()?;
    
    Ok(())
}