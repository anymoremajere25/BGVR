use std::collections::HashMap;
use std::path::PathBuf;
use std::fs::File;
use std::io::Write;
use anyhow::{Context, Result};
use clap::Parser;
use log::{info, warn, error};
use rayon::prelude::*;
use rust_htslib::bam::{self, Read, Record};
use serde::{Deserialize, Serialize};
use statrs::distribution::{Poisson, Continuous};
use itertools::Itertools;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input BAM file path
    #[arg(short, long)]
    input: PathBuf,
    
    /// Output peaks file path
    #[arg(short, long, default_value = "peaks.json")]
    output: PathBuf,
    
    /// Window size for peak calling (bp)
    #[arg(short, long, default_value_t = 200)]
    window_size: u64,
    
    /// Minimum coverage threshold
    #[arg(short, long, default_value_t = 5.0)]
    min_coverage: f64,
    
    /// P-value threshold for peak significance
    #[arg(short, long, default_value_t = 0.05)]
    pvalue_threshold: f64,
    
    /// Number of threads to use
    #[arg(short, long, default_value_t = 4)]
    threads: usize,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
struct Peak {
    chrom: String,
    start: u64,
    end: u64,
    summit: u64,
    coverage: f64,
    pvalue: f64,
    qvalue: f64,
    fold_enrichment: f64,
}

#[derive(Debug)]
struct GenomicInterval {
    chrom: String,
    start: u64,
    end: u64,
    coverage: f64,
}

struct PeakCaller {
    window_size: u64,
    min_coverage: f64,
    pvalue_threshold: f64,
}

impl PeakCaller {
    fn new(window_size: u64, min_coverage: f64, pvalue_threshold: f64) -> Self {
        Self {
            window_size,
            min_coverage,
            pvalue_threshold,
        }
    }

    fn calculate_coverage(&self, reads: &[(String, u64, u64)]) -> HashMap<String, Vec<(u64, f64)>> {
        let mut coverage_map: HashMap<String, HashMap<u64, u64>> = HashMap::new();
        
        // Count reads per window
        for (chrom, start, end) in reads {
            let chrom_map = coverage_map.entry(chrom.clone()).or_insert_with(HashMap::new);
            
            let window_start = (start / self.window_size) * self.window_size;
            let window_end = (end / self.window_size) * self.window_size;
            
            for window in (window_start..=window_end).step_by(self.window_size as usize) {
                *chrom_map.entry(window).or_insert(0) += 1;
            }
        }
        
        // Convert to sorted vectors
        coverage_map
            .into_iter()
            .map(|(chrom, windows)| {
                let mut sorted_windows: Vec<(u64, f64)> = windows
                    .into_iter()
                    .map(|(pos, count)| (pos, count as f64))
                    .collect();
                sorted_windows.sort_by_key(|&(pos, _)| pos);
                (chrom, sorted_windows)
            })
            .collect()
    }

    fn call_peaks(&self, coverage_data: HashMap<String, Vec<(u64, f64)>>) -> Result<Vec<Peak>> {
        let mut all_peaks = Vec::new();
        
        for (chrom, windows) in coverage_data {
            if windows.is_empty() {
                continue;
            }
            
            // Calculate background rate (mean coverage)
            let total_coverage: f64 = windows.iter().map(|(_, cov)| cov).sum();
            let background_rate = total_coverage / windows.len() as f64;
            
            if background_rate < 0.1 {
                warn!("Very low background coverage for chromosome {}: {:.2}", chrom, background_rate);
                continue;
            }
            
            info!("Processing chromosome {} with background rate: {:.2}", chrom, background_rate);
            
            // Find peaks using Poisson distribution
            let poisson = Poisson::new(background_rate)
                .context("Failed to create Poisson distribution")?;
            
            let mut chromosome_peaks = Vec::new();
            
            for (pos, coverage) in &windows {
                if *coverage >= self.min_coverage {
                    // Calculate p-value using Poisson distribution
                    let pvalue = 1.0 - poisson.cdf(*coverage - 1.0);
                    
                    if pvalue <= self.pvalue_threshold {
                        let fold_enrichment = coverage / background_rate.max(1.0);
                        
                        chromosome_peaks.push(Peak {
                            chrom: chrom.clone(),
                            start: *pos,
                            end: pos + self.window_size,
                            summit: pos + self.window_size / 2,
                            coverage: *coverage,
                            pvalue,
                            qvalue: pvalue, // Simplified - in practice, use Benjamini-Hochberg
                            fold_enrichment,
                        });
                    }
                }
            }
            
            // Merge nearby peaks
            chromosome_peaks = self.merge_nearby_peaks(chromosome_peaks);
            all_peaks.extend(chromosome_peaks);
        }
        
        // Sort peaks by significance
        all_peaks.sort_by(|a, b| a.pvalue.partial_cmp(&b.pvalue).unwrap_or(std::cmp::Ordering::Equal));
        
        // Apply multiple testing correction (simplified Benjamini-Hochberg)
        let n_tests = all_peaks.len() as f64;
        for (i, peak) in all_peaks.iter_mut().enumerate() {
            let rank = (i + 1) as f64;
            peak.qvalue = (peak.pvalue * n_tests / rank).min(1.0);
        }
        
        Ok(all_peaks)
    }

    fn merge_nearby_peaks(&self, mut peaks: Vec<Peak>) -> Vec<Peak> {
        if peaks.is_empty() {
            return peaks;
        }
        
        peaks.sort_by_key(|p| p.start);
        let mut merged = Vec::new();
        let mut current = peaks[0].clone();
        
        for peak in peaks.into_iter().skip(1) {
            if peak.start <= current.end + self.window_size {
                // Merge peaks - keep the one with higher coverage as summit
                if peak.coverage > current.coverage {
                    current.summit = peak.summit;
                    current.coverage = peak.coverage;
                    current.pvalue = peak.pvalue.min(current.pvalue);
                }
                current.end = peak.end.max(current.end);
                current.fold_enrichment = current.fold_enrichment.max(peak.fold_enrichment);
            } else {
                merged.push(current);
                current = peak;
            }
        }
        merged.push(current);
        merged
    }
}

fn read_bam_file(bam_path: &PathBuf) -> Result<Vec<(String, u64, u64)>> {
    info!("Reading BAM file: {:?}", bam_path);
    
    let mut bam_reader = bam::Reader::from_path(bam_path)
        .with_context(|| format!("Failed to open BAM file: {:?}", bam_path))?;
    
    let header = bam_reader.header().clone();
    let mut reads = Vec::new();
    let mut record_count = 0;
    let mut valid_records = 0;

    for result in bam_reader.records() {
        record_count += 1;
        
        match result {
            Ok(record) => {
                if record.is_unmapped() || record.is_secondary() || record.is_duplicate() {
                    continue;
                }
                
                let tid = record.tid();
                if tid < 0 {
                    continue;
                }
                
                let chrom_name = match header.tid2name(tid as u32) {
                    Some(name) => String::from_utf8_lossy(name).to_string(),
                    None => {
                        warn!("Unknown chromosome ID: {}", tid);
                        continue;
                    }
                };
                
                let start = record.pos() as u64;
                let end = record.cigar_cached()
                    .map(|cigar| start + cigar.end_pos() as u64)
                    .unwrap_or(start + record.seq_len() as u64);
                
                reads.push((chrom_name, start, end));
                valid_records += 1;
            }
            Err(e) => {
                error!("Error reading BAM record {}: {}", record_count, e);
            }
        }
        
        if record_count % 100000 == 0 {
            info!("Processed {} records, {} valid", record_count, valid_records);
        }
    }
    
    info!("Finished reading BAM: {} total records, {} valid alignments", record_count, valid_records);
    Ok(reads)
}

fn write_peaks(peaks: &[Peak], output_path: &PathBuf) -> Result<()> {
    info!("Writing {} peaks to {:?}", peaks.len(), output_path);
    
    let output = serde_json::to_string_pretty(peaks)
        .context("Failed to serialize peaks to JSON")?;
    
    let mut file = File::create(output_path)
        .with_context(|| format!("Failed to create output file: {:?}", output_path))?;
    
    file.write_all(output.as_bytes())
        .with_context(|| format!("Failed to write to output file: {:?}", output_path))?;
    
    Ok(())
}

fn main() -> Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    
    let args = Args::parse();
    
    info!("Starting peak calling with parameters:");
    info!("  Input: {:?}", args.input);
    info!("  Output: {:?}", args.output);
    info!("  Window size: {} bp", args.window_size);
    info!("  Min coverage: {}", args.min_coverage);
    info!("  P-value threshold: {}", args.pvalue_threshold);
    info!("  Threads: {}", args.threads);
    
    // Set number of threads for rayon
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .context("Failed to initialize thread pool")?;
    
    // Validate input file exists
    if !args.input.exists() {
        anyhow::bail!("Input BAM file does not exist: {:?}", args.input);
    }
    
    // Read BAM file
    let reads = read_bam_file(&args.input)?;
    
    if reads.is_empty() {
        anyhow::bail!("No valid reads found in BAM file");
    }
    
    // Initialize peak caller
    let peak_caller = PeakCaller::new(
        args.window_size,
        args.min_coverage,
        args.pvalue_threshold,
    );
    
    // Calculate coverage
    info!("Calculating coverage across genome...");
    let coverage_data = peak_caller.calculate_coverage(&reads);
    
    // Call peaks
    info!("Calling peaks...");
    let peaks = peak_caller.call_peaks(coverage_data)?;
    
    info!("Found {} significant peaks", peaks.len());
    
    // Write results
    write_peaks(&peaks, &args.output)?;
    
    // Print summary statistics
    if !peaks.is_empty() {
        let avg_coverage: f64 = peaks.iter().map(|p| p.coverage).sum::<f64>() / peaks.len() as f64;
        let min_pvalue = peaks.iter().map(|p| p.pvalue).fold(1.0, f64::min);
        let max_fold_enrichment = peaks.iter().map(|p| p.fold_enrichment).fold(0.0, f64::max);
        
        info!("Peak calling summary:");
        info!("  Total peaks: {}", peaks.len());
        info!("  Average coverage: {:.2}", avg_coverage);
        info!("  Best p-value: {:.2e}", min_pvalue);
        info!("  Max fold enrichment: {:.2}", max_fold_enrichment);
    }
    
    info!("Peak calling completed successfully!");
    Ok(())
}