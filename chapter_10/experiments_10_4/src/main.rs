use std::collections::HashMap;
use std::path::PathBuf;
use std::fs::File;
use std::io::Write;
use std::sync::Arc;
use anyhow::{Context, Result, bail};
use clap::Parser;
use log::{info, warn, error, debug};
use rayon::prelude::*;
use rust_htslib::bam::{self, Read, Record, HeaderView};
use serde::{Deserialize, Serialize};
use statrs::statistics::{Statistics, OrderStatistics};
use indicatif::{ProgressBar, ProgressStyle};
use chrono::{DateTime, Utc};
use ndarray::{Array1, Axis};

#[cfg(feature = "fft-acceleration")]
use rustfft::{FftPlanner, num_complex::Complex};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input BAM file path
    #[arg(short, long)]
    input: PathBuf,
    
    /// Output shift estimate file (JSON format)
    #[arg(short, long, default_value = "shift_estimate.json")]
    output: PathBuf,
    
    /// Output detailed correlation file
    #[arg(long)]
    correlation_output: Option<PathBuf>,
    
    /// Output coverage profiles
    #[arg(long)]
    coverage_output: Option<PathBuf>,
    
    /// Maximum shift range to search (bp)
    #[arg(long, default_value_t = 500)]
    max_shift: u32,
    
    /// Minimum mapping quality
    #[arg(long, default_value_t = 10)]
    min_mapq: u8,
    
    /// Chromosome to analyze (default: all)
    #[arg(long)]
    chromosome: Option<String>,
    
    /// Start position for analysis (1-based)
    #[arg(long)]
    start_pos: Option<u64>,
    
    /// End position for analysis (1-based)
    #[arg(long)]
    end_pos: Option<u64>,
    
    /// Bin size for coverage calculation (bp)
    #[arg(long, default_value_t = 1)]
    bin_size: u32,
    
    /// Sampling factor (analyze every Nth read)
    #[arg(long, default_value_t = 1)]
    sampling_factor: u32,
    
    /// Number of threads to use
    #[arg(short, long, default_value_t = 4)]
    threads: usize,
    
    /// Enable verbose logging
    #[arg(short, long)]
    verbose: bool,
    
    /// Use FFT acceleration for cross-correlation
    #[arg(long)]
    use_fft: bool,
    
    /// Smoothing window size for coverage profiles
    #[arg(long, default_value_t = 5)]
    smoothing_window: usize,
}

#[derive(Debug, Serialize, Deserialize)]
struct ShiftEstimate {
    estimated_shift: i64,
    confidence_score: f64,
    correlation_peak: f64,
    background_correlation: f64,
    signal_to_noise_ratio: f64,
    phantom_peak_shift: Option<i64>,
    read_length_estimate: f64,
    analysis_region: AnalysisRegion,
    quality_metrics: QualityMetrics,
    processing_stats: ProcessingStats,
}

#[derive(Debug, Serialize, Deserialize)]
struct AnalysisRegion {
    chromosome: Option<String>,
    start: Option<u64>,
    end: Option<u64>,
    total_length: u64,
}

#[derive(Debug, Serialize, Deserialize)]
struct QualityMetrics {
    normalized_strand_correlation: f64,
    relative_strand_correlation: f64,
    reads_processed: u64,
    forward_reads: u64,
    reverse_reads: u64,
    correlation_at_read_length: f64,
}

#[derive(Debug, Serialize, Deserialize)]
struct ProcessingStats {
    analysis_time_seconds: f64,
    memory_usage_mb: f64,
    bin_size: u32,
    max_shift_searched: u32,
    fft_acceleration_used: bool,
}

#[derive(Debug, Serialize, Deserialize)]
struct CorrelationProfile {
    shift_values: Vec<i64>,
    correlation_values: Vec<f64>,
    smoothed_correlation: Vec<f64>,
}

struct CoverageProfile {
    forward_coverage: Array1<f64>,
    reverse_coverage: Array1<f64>,
    positions: Vec<u64>,
    bin_size: u32,
}

struct CrossCorrelationAnalyzer {
    args: Arc<Args>,
}

impl CrossCorrelationAnalyzer {
    fn new(args: Args) -> Self {
        Self {
            args: Arc::new(args),
        }
    }
    
    fn calculate_coverage_profiles(&self, bam_path: &PathBuf) -> Result<CoverageProfile> {
        info!("Calculating strand-specific coverage profiles...");
        
        let mut bam_reader = bam::Reader::from_path(bam_path)
            .with_context(|| format!("Failed to open BAM file: {:?}", bam_path))?;
        
        let header = bam_reader.header().clone();
        
        // Determine analysis region
        let (chrom_name, start_pos, end_pos, total_length) = 
            self.determine_analysis_region(&header)?;
        
        info!("Analysis region: {}:{}-{} ({} bp)", 
              chrom_name.as_ref().unwrap_or(&"all".to_string()), 
              start_pos.unwrap_or(1), 
              end_pos.unwrap_or(total_length), 
              total_length);
        
        // Calculate number of bins
        let num_bins = (total_length as f64 / self.args.bin_size as f64).ceil() as usize;
        let mut forward_coverage = Array1::<f64>::zeros(num_bins);
        let mut reverse_coverage = Array1::<f64>::zeros(num_bins);
        let mut positions = Vec::with_capacity(num_bins);
        
        // Generate position vector
        for i in 0..num_bins {
            let pos = start_pos.unwrap_or(0) + (i as u64 * self.args.bin_size as u64);
            positions.push(pos);
        }
        
        // Set up progress tracking
        let pb = ProgressBar::new_spinner();
        pb.set_style(ProgressStyle::default_spinner()
            .template("{spinner:.green} Processing reads: {pos}").unwrap());
        
        let mut reads_processed = 0u64;
        let mut forward_reads = 0u64;
        let mut reverse_reads = 0u64;
        
        // Process reads
        for result in bam_reader.records() {
            let record = result?;
            reads_processed += 1;
            
            if reads_processed % 10000 == 0 {
                pb.set_position(reads_processed);
            }
            
            // Apply sampling
            if reads_processed % self.args.sampling_factor as u64 != 0 {
                continue;
            }
            
            // Filter reads
            if !self.should_include_read(&record, &chrom_name, start_pos, end_pos) {
                continue;
            }
            
            // Calculate bin index
            let read_pos = record.pos() as u64;
            let relative_pos = read_pos - start_pos.unwrap_or(0);
            let bin_index = (relative_pos / self.args.bin_size as u64) as usize;
            
            if bin_index < num_bins {
                if record.is_reverse() {
                    reverse_coverage[bin_index] += 1.0;
                    reverse_reads += 1;
                } else {
                    forward_coverage[bin_index] += 1.0;
                    forward_reads += 1;
                }
            }
        }
        
        pb.finish_with_message(format!("Processed {} reads ({} forward, {} reverse)", 
                                      reads_processed, forward_reads, reverse_reads));
        
        info!("Coverage calculation complete");
        debug!("Forward reads: {}, Reverse reads: {}", forward_reads, reverse_reads);
        
        // Apply smoothing if requested
        if self.args.smoothing_window > 1 {
            info!("Applying smoothing with window size {}", self.args.smoothing_window);
            forward_coverage = self.smooth_coverage(&forward_coverage);
            reverse_coverage = self.smooth_coverage(&reverse_coverage);
        }
        
        Ok(CoverageProfile {
            forward_coverage,
            reverse_coverage,
            positions,
            bin_size: self.args.bin_size,
        })
    }
    
    fn determine_analysis_region(&self, header: &HeaderView) -> Result<(Option<String>, Option<u64>, Option<u64>, u64)> {
        if let Some(ref chrom) = self.args.chromosome {
            // Analyze specific chromosome
            let tid = header.tid(chrom.as_bytes())
                .with_context(|| format!("Chromosome '{}' not found in BAM header", chrom))?;
            
            let chrom_length = header.target_len(tid)
                .with_context(|| format!("Failed to get length for chromosome '{}'", chrom))? as u64;
            
            let start = self.args.start_pos.unwrap_or(1);
            let end = self.args.end_pos.unwrap_or(chrom_length);
            
            if start >= end {
                bail!("Start position {} must be less than end position {}", start, end);
            }
            
            let length = end - start + 1;
            Ok((Some(chrom.clone()), Some(start), Some(end), length))
        } else {
            // Analyze entire genome
            let total_length: u64 = (0..header.target_count())
                .map(|i| header.target_len(i).unwrap_or(0) as u64)
                .sum();
            
            Ok((None, None, None, total_length))
        }
    }
    
    fn should_include_read(&self, record: &Record, target_chrom: &Option<String>, 
                          start_pos: Option<u64>, end_pos: Option<u64>) -> bool {
        // Basic quality filters
        if record.is_unmapped() || 
           record.is_secondary() || 
           record.is_duplicate() ||
           record.mapq() < self.args.min_mapq {
            return false;
        }
        
        // Chromosome filter
        if let Some(ref target) = target_chrom {
            let tid = record.tid();
            if tid < 0 {
                return false;
            }
            // Note: This is a simplified check - in practice, we'd verify chromosome name
        }
        
        // Position filter
        if let (Some(start), Some(end)) = (start_pos, end_pos) {
            let read_pos = record.pos() as u64;
            if read_pos < start || read_pos > end {
                return false;
            }
        }
        
        true
    }
    
    fn smooth_coverage(&self, coverage: &Array1<f64>) -> Array1<f64> {
        let window = self.args.smoothing_window;
        let mut smoothed = Array1::<f64>::zeros(coverage.len());
        
        for i in 0..coverage.len() {
            let start = i.saturating_sub(window / 2);
            let end = (i + window / 2 + 1).min(coverage.len());
            
            let sum: f64 = coverage.slice(ndarray::s![start..end]).sum();
            let count = end - start;
            smoothed[i] = sum / count as f64;
        }
        
        smoothed
    }
    
    #[cfg(feature = "fft-acceleration")]
    fn cross_correlate_fft(&self, forward: &Array1<f64>, reverse: &Array1<f64>) -> Result<Array1<f64>> {
        info!("Computing cross-correlation using FFT acceleration...");
        
        let n = forward.len();
        let fft_size = (2 * n).next_power_of_two();
        
        let mut planner = FftPlanner::<f64>::new();
        let fft = planner.plan_fft_forward(fft_size);
        let ifft = planner.plan_fft_inverse(fft_size);
        
        // Prepare forward signal
        let mut forward_complex: Vec<Complex<f64>> = forward
            .iter()
            .map(|&x| Complex::new(x, 0.0))
            .chain(std::iter::repeat(Complex::new(0.0, 0.0)))
            .take(fft_size)
            .collect();
        
        // Prepare reverse signal (conjugated for correlation)
        let mut reverse_complex: Vec<Complex<f64>> = reverse
            .iter()
            .map(|&x| Complex::new(x, 0.0))
            .chain(std::iter::repeat(Complex::new(0.0, 0.0)))
            .take(fft_size)
            .collect();
        
        // Reverse and conjugate for correlation
        reverse_complex.reverse();
        for c in &mut reverse_complex {
            *c = c.conj();
        }
        
        // Forward FFT
        fft.process(&mut forward_complex);
        fft.process(&mut reverse_complex);
        
        // Element-wise multiplication
        for i in 0..fft_size {
            forward_complex[i] *= reverse_complex[i];
        }
        
        // Inverse FFT
        ifft.process(&mut forward_complex);
        
        // Extract real parts and normalize
        let correlation: Vec<f64> = forward_complex
            .iter()
            .take(2 * n - 1)
            .map(|c| c.re / fft_size as f64)
            .collect();
        
        Ok(Array1::from_vec(correlation))
    }
    
    fn cross_correlate_direct(&self, forward: &Array1<f64>, reverse: &Array1<f64>) -> Result<Array1<f64>> {
        info!("Computing cross-correlation using direct method...");
        
        let n = forward.len();
        let max_lag = (self.args.max_shift as usize / self.args.bin_size as usize).min(n / 2);
        
        let pb = ProgressBar::new(2 * max_lag as u64 + 1);
        pb.set_style(ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos:>7}/{len:7} {msg}")
            .unwrap());
        
        let correlations: Vec<f64> = (-max_lag as i32..=max_lag as i32)
            .into_par_iter()
            .map(|lag| {
                pb.inc(1);
                
                let mut correlation = 0.0;
                let mut count = 0;
                
                for i in 0..n {
                    let j = i as i32 + lag;
                    if j >= 0 && (j as usize) < n {
                        correlation += forward[i] * reverse[j as usize];
                        count += 1;
                    }
                }
                
                if count > 0 {
                    correlation / count as f64
                } else {
                    0.0
                }
            })
            .collect();
        
        pb.finish_with_message("Cross-correlation completed");
        
        Ok(Array1::from_vec(correlations))
    }
    
    fn compute_cross_correlation(&self, coverage: &CoverageProfile) -> Result<CorrelationProfile> {
        info!("Computing cross-correlation between forward and reverse strands...");
        
        // Choose correlation method
        let correlation_values = if self.args.use_fft && cfg!(feature = "fft-acceleration") {
            #[cfg(feature = "fft-acceleration")]
            {
                self.cross_correlate_fft(&coverage.forward_coverage, &coverage.reverse_coverage)?
            }
            #[cfg(not(feature = "fft-acceleration"))]
            {
                warn!("FFT acceleration not available, falling back to direct method");
                self.cross_correlate_direct(&coverage.forward_coverage, &coverage.reverse_coverage)?
            }
        } else {
            self.cross_correlate_direct(&coverage.forward_coverage, &coverage.reverse_coverage)?
        };
        
        // Generate shift values
        let max_lag = (self.args.max_shift as usize / self.args.bin_size as usize)
            .min(coverage.forward_coverage.len() / 2);
        
        let shift_values: Vec<i64> = (-max_lag as i64..=max_lag as i64)
            .map(|lag| lag * self.args.bin_size as i64)
            .collect();
        
        // Apply additional smoothing to correlation profile
        let smoothed_correlation = self.smooth_correlation_profile(&correlation_values);
        
        Ok(CorrelationProfile {
            shift_values,
            correlation_values: correlation_values.to_vec(),
            smoothed_correlation,
        })
    }
    
    fn smooth_correlation_profile(&self, correlation: &Array1<f64>) -> Vec<f64> {
        let window = 3; // Small smoothing window for correlation
        let mut smoothed = vec![0.0; correlation.len()];
        
        for i in 0..correlation.len() {
            let start = i.saturating_sub(window / 2);
            let end = (i + window / 2 + 1).min(correlation.len());
            
            let sum: f64 = correlation.slice(ndarray::s![start..end]).sum();
            let count = end - start;
            smoothed[i] = sum / count as f64;
        }
        
        smoothed
    }
    
    fn analyze_correlation_profile(&self, profile: &CorrelationProfile) -> Result<ShiftEstimate> {
        info!("Analyzing correlation profile to estimate fragment shift...");
        
        let correlations = &profile.smoothed_correlation;
        
        // Find the peak correlation (maximum)
        let (peak_index, &peak_correlation) = correlations
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .context("Failed to find correlation peak")?;
        
        let estimated_shift = profile.shift_values[peak_index];
        
        // Calculate background correlation (mean of correlation values)
        let background_correlation = correlations.iter().sum::<f64>() / correlations.len() as f64;
        
        // Calculate signal-to-noise ratio
        let correlation_std = {
            let mean = background_correlation;
            let variance = correlations.iter()
                .map(|&x| (x - mean).powi(2))
                .sum::<f64>() / correlations.len() as f64;
            variance.sqrt()
        };
        
        let signal_to_noise_ratio = if correlation_std > 0.0 {
            (peak_correlation - background_correlation) / correlation_std
        } else {
            0.0
        };
        
        // Find phantom peak (read length artifact)
        let phantom_peak_shift = self.find_phantom_peak(profile);
        
        // Calculate quality metrics
        let quality_metrics = self.calculate_quality_metrics(profile, estimated_shift)?;
        
        // Calculate confidence score based on peak prominence
        let confidence_score = self.calculate_confidence_score(
            peak_correlation, 
            background_correlation, 
            signal_to_noise_ratio
        );
        
        info!("Fragment shift estimate: {} bp", estimated_shift);
        info!("Peak correlation: {:.4}", peak_correlation);
        info!("Background correlation: {:.4}", background_correlation);
        info!("Signal-to-noise ratio: {:.2}", signal_to_noise_ratio);
        info!("Confidence score: {:.3}", confidence_score);
        
        Ok(ShiftEstimate {
            estimated_shift,
            confidence_score,
            correlation_peak: peak_correlation,
            background_correlation,
            signal_to_noise_ratio,
            phantom_peak_shift,
            read_length_estimate: 75.0, // Placeholder - would estimate from data
            analysis_region: AnalysisRegion {
                chromosome: self.args.chromosome.clone(),
                start: self.args.start_pos,
                end: self.args.end_pos,
                total_length: 0, // Would be filled from coverage analysis
            },
            quality_metrics,
            processing_stats: ProcessingStats {
                analysis_time_seconds: 0.0, // Would be filled by caller
                memory_usage_mb: 0.0, // Would be estimated
                bin_size: self.args.bin_size,
                max_shift_searched: self.args.max_shift,
                fft_acceleration_used: self.args.use_fft && cfg!(feature = "fft-acceleration"),
            },
        })
    }
    
    fn find_phantom_peak(&self, profile: &CorrelationProfile) -> Option<i64> {
        // Look for secondary peak around read length (typically 36-150 bp)
        let read_length_range = 36..=150;
        let mut best_phantom = None;
        let mut best_correlation = 0.0;
        
        for (i, &shift) in profile.shift_values.iter().enumerate() {
            if read_length_range.contains(&shift.abs()) {
                let correlation = profile.smoothed_correlation[i];
                if correlation > best_correlation {
                    best_correlation = correlation;
                    best_phantom = Some(shift);
                }
            }
        }
        
        best_phantom
    }
    
    fn calculate_quality_metrics(&self, profile: &CorrelationProfile, estimated_shift: i64) -> Result<QualityMetrics> {
        // Simplified quality metrics calculation
        // In practice, these would be more sophisticated (e.g., NSC, RSC from phantompeakqualtools)
        
        let correlations = &profile.smoothed_correlation;
        let max_correlation = correlations.iter().fold(0.0, |a, &b| a.max(b));
        let min_correlation = correlations.iter().fold(f64::MAX, |a, &b| a.min(b));
        
        // Normalized Strand Correlation (NSC) approximation
        let nsc = if min_correlation != 0.0 {
            max_correlation / min_correlation
        } else {
            max_correlation
        };
        
        // Relative Strand Correlation (RSC) approximation
        let background = correlations.iter().sum::<f64>() / correlations.len() as f64;
        let rsc = if background != 0.0 {
            (max_correlation - background) / (background - min_correlation).max(0.001)
        } else {
            0.0
        };
        
        // Find correlation at read length
        let read_length_correlation = profile.shift_values
            .iter()
            .position(|&shift| shift.abs() == 75) // Assume 75bp read length
            .map(|i| correlations[i])
            .unwrap_or(0.0);
        
        Ok(QualityMetrics {
            normalized_strand_correlation: nsc,
            relative_strand_correlation: rsc,
            reads_processed: 0, // Would be filled from coverage analysis
            forward_reads: 0,   // Would be filled from coverage analysis
            reverse_reads: 0,   // Would be filled from coverage analysis
            correlation_at_read_length: read_length_correlation,
        })
    }
    
    fn calculate_confidence_score(&self, peak_correlation: f64, background_correlation: f64, snr: f64) -> f64 {
        // Combine multiple factors into a confidence score (0-1 scale)
        let peak_strength = (peak_correlation - background_correlation) / peak_correlation.max(0.001);
        let snr_normalized = (snr / 10.0).min(1.0).max(0.0); // Normalize SNR to 0-1
        
        // Weighted combination
        0.6 * peak_strength + 0.4 * snr_normalized
    }
    
    fn write_outputs(&self, shift_estimate: &ShiftEstimate, 
                    correlation_profile: &Option<CorrelationProfile>) -> Result<()> {
        // Write main shift estimate
        info!("Writing shift estimate to {:?}", self.args.output);
        let output = serde_json::to_string_pretty(shift_estimate)
            .context("Failed to serialize shift estimate")?;
        
        let mut file = File::create(&self.args.output)
            .with_context(|| format!("Failed to create output file: {:?}", self.args.output))?;
        
        file.write_all(output.as_bytes())
            .with_context(|| format!("Failed to write to output file: {:?}", self.args.output))?;
        
        // Write correlation profile if requested
        if let (Some(output_path), Some(profile)) = (&self.args.correlation_output, correlation_profile) {
            info!("Writing correlation profile to {:?}", output_path);
            let correlation_output = serde_json::to_string_pretty(profile)
                .context("Failed to serialize correlation profile")?;
            
            let mut file = File::create(output_path)
                .with_context(|| format!("Failed to create correlation output file: {:?}", output_path))?;
            
            file.write_all(correlation_output.as_bytes())
                .with_context(|| format!("Failed to write correlation file: {:?}", output_path))?;
        }
        
        Ok(())
    }
    
    pub fn run(&mut self) -> Result<()> {
        let start_time = Utc::now();
        info!("Starting cross-correlation fragment shift analysis at {}", 
              start_time.format("%Y-%m-%d %H:%M:%S"));
        
        // Validate input file
        if !self.args.input.exists() {
            bail!("Input BAM file does not exist: {:?}", self.args.input);
        }
        
        // Calculate coverage profiles
        let coverage = self.calculate_coverage_profiles(&self.args.input)?;
        
        // Compute cross-correlation
        let correlation_profile = self.compute_cross_correlation(&coverage)?;
        
        // Analyze correlation to estimate shift
        let mut shift_estimate = self.analyze_correlation_profile(&correlation_profile)?;
        
        // Fill in processing statistics
        let end_time = Utc::now();
        let duration = end_time.signed_duration_since(start_time);
        shift_estimate.processing_stats.analysis_time_seconds = 
            duration.num_milliseconds() as f64 / 1000.0;
        
        // Write outputs
        self.write_outputs(&shift_estimate, &Some(correlation_profile))?;
        
        info!("Fragment shift analysis completed in {:.1} seconds", 
              shift_estimate.processing_stats.analysis_time_seconds);
        
        // Print summary
        println!("\n=== Fragment Shift Analysis Summary ===");
        println!("Estimated fragment shift: {} bp", shift_estimate.estimated_shift);
        println!("Confidence score: {:.3}", shift_estimate.confidence_score);
        println!("Signal-to-noise ratio: {:.2}", shift_estimate.signal_to_noise_ratio);
        println!("Peak correlation: {:.4}", shift_estimate.correlation_peak);
        if let Some(phantom_shift) = shift_estimate.phantom_peak_shift {
            println!("Phantom peak shift: {} bp", phantom_shift);
        }
        println!("Analysis time: {:.1} seconds", shift_estimate.processing_stats.analysis_time_seconds);
        println!("======================================");
        
        Ok(())
    }
}

fn main() -> Result<()> {
    let args = Args::parse();
    
    // Initialize logging
    let log_level = if args.verbose { "debug" } else { "info" };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(log_level)).init();
    
    info!("Fragment shift estimator starting...");
    info!("Parameters:");
    info!("  Input BAM: {:?}", args.input);
    info!("  Output: {:?}", args.output);
    info!("  Max shift: {} bp", args.max_shift);
    info!("  Min MAPQ: {}", args.min_mapq);
    info!("  Bin size: {} bp", args.bin_size);
    info!("  Chromosome: {:?}", args.chromosome.as_ref().unwrap_or(&"all".to_string()));
    info!("  Use FFT: {}", args.use_fft);
    info!("  Threads: {}", args.threads);
    
    // Set up thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .context("Failed to initialize thread pool")?;
    
    // Run analysis
    let mut analyzer = CrossCorrelationAnalyzer::new(args);
    analyzer.run()?;
    
    Ok(())
}