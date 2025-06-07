use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

// Optional imports for enhanced functionality
#[cfg(feature = "parallel")]
use rayon::prelude::*;

#[cfg(feature = "compression")]
use flate2::read::GzDecoder;

#[cfg(feature = "json")]
use serde::{Deserialize, Serialize};

#[cfg(feature = "bio-utils")]
use bio::io::bed;

/// Represents a genomic interval with coverage data
#[derive(Debug, Clone)]
#[cfg_attr(feature = "json", derive(Serialize, Deserialize))]
pub struct CoverageInterval {
    pub chromosome: String,
    pub start: u32,
    pub end: u32,
    pub coverage: f64,
}

/// Represents a called peak
#[derive(Debug, Clone)]
#[cfg_attr(feature = "json", derive(Serialize, Deserialize))]
pub struct Peak {
    pub chromosome: String,
    pub start: u32,
    pub end: u32,
    pub signal: f64,
    pub background: f64,
    pub fold_enrichment: f64,
}

/// Configuration parameters for peak calling
#[derive(Debug, Clone)]
pub struct PeakCallingConfig {
    /// Minimum signal-to-background ratio to call a peak
    pub min_fold_enrichment: f64,
    /// Window size for local background calculation (bp)
    pub background_window: u32,
    /// Minimum peak length (bp)
    pub min_peak_length: u32,
    /// Merge peaks within this distance (bp)
    pub merge_distance: u32,
}

impl Default for PeakCallingConfig {
    fn default() -> Self {
        Self {
            min_fold_enrichment: 2.0,
            background_window: 10000,
            min_peak_length: 150,
            merge_distance: 200,
        }
    }
}

/// Main peak calling module
pub struct PeakCaller {
    /// Coverage data organized by chromosome
    coverage_data: BTreeMap<String, Vec<CoverageInterval>>,
    config: PeakCallingConfig,
}

impl PeakCaller {
    /// Create a new peak caller with default configuration
    pub fn new() -> Self {
        Self {
            coverage_data: BTreeMap::new(),
            config: PeakCallingConfig::default(),
        }
    }

    /// Create a new peak caller with custom configuration
    pub fn with_config(config: PeakCallingConfig) -> Self {
        Self {
            coverage_data: BTreeMap::new(),
            config,
        }
    }

    /// Load coverage data from a BedGraph file (supports .gz compression)
    /// Format: chromosome\tstart\tend\tcoverage
    pub fn load_coverage_from_bedgraph<P: AsRef<Path>>(&mut self, file_path: P) -> Result<(), Box<dyn std::error::Error>> {
        let file = File::open(&file_path)?;
        
        // Handle compressed files if compression feature is enabled
        #[cfg(feature = "compression")]
        let reader: Box<dyn BufRead> = if file_path.as_ref().extension()
            .and_then(|s| s.to_str())
            .map(|s| s.to_lowercase())
            .as_deref() == Some("gz") {
            let gz_decoder = GzDecoder::new(file);
            Box::new(BufReader::new(gz_decoder))
        } else {
            Box::new(BufReader::new(file))
        };
        
        #[cfg(not(feature = "compression"))]
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let line = line?;
            if line.starts_with('#') || line.starts_with("track") || line.trim().is_empty() {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 4 {
                continue;
            }

            let chromosome = fields[0].to_string();
            let start: u32 = fields[1].parse()?;
            let end: u32 = fields[2].parse()?;
            let coverage: f64 = fields[3].parse()?;

            let interval = CoverageInterval {
                chromosome: chromosome.clone(),
                start,
                end,
                coverage,
            };

            self.coverage_data
                .entry(chromosome)
                .or_insert_with(Vec::new)
                .push(interval);
        }

        // Sort intervals by start position for each chromosome
        for intervals in self.coverage_data.values_mut() {
            intervals.sort_by_key(|interval| interval.start);
        }

        Ok(())
    }

    /// Add a single coverage interval
    pub fn add_coverage_interval(&mut self, interval: CoverageInterval) {
        let chromosome = interval.chromosome.clone();
        self.coverage_data
            .entry(chromosome)
            .or_insert_with(Vec::new)
            .push(interval);
    }

    /// Calculate local background coverage for a given position
    fn calculate_local_background(&self, chromosome: &str, position: u32) -> f64 {
        let Some(intervals) = self.coverage_data.get(chromosome) else {
            return 0.0;
        };

        let window_start = position.saturating_sub(self.config.background_window / 2);
        let window_end = position + self.config.background_window / 2;

        let mut total_coverage = 0.0;
        let mut total_length = 0u32;

        for interval in intervals {
            // Skip intervals that don't overlap with our background window
            if interval.end <= window_start || interval.start >= window_end {
                continue;
            }

            // Calculate overlap between interval and background window
            let overlap_start = interval.start.max(window_start);
            let overlap_end = interval.end.min(window_end);
            let overlap_length = overlap_end - overlap_start;

            total_coverage += interval.coverage * overlap_length as f64;
            total_length += overlap_length;
        }

        if total_length > 0 {
            total_coverage / total_length as f64
        } else {
            0.0
        }
    }

    /// Get coverage at a specific position using binary search
    fn get_coverage_at_position(&self, chromosome: &str, position: u32) -> f64 {
        let Some(intervals) = self.coverage_data.get(chromosome) else {
            return 0.0;
        };

        // Binary search for the interval containing this position
        match intervals.binary_search_by(|interval| {
            if position < interval.start {
                std::cmp::Ordering::Greater
            } else if position >= interval.end {
                std::cmp::Ordering::Less
            } else {
                std::cmp::Ordering::Equal
            }
        }) {
            Ok(index) => intervals[index].coverage,
            Err(_) => 0.0,
        }
    }

    /// Call peaks for all chromosomes (with optional parallelization)
    pub fn call_peaks(&self) -> Vec<Peak> {
        #[cfg(feature = "parallel")]
        let chromosome_peaks: Vec<Vec<Peak>> = self.coverage_data
            .par_iter()
            .map(|(chromosome, intervals)| {
                self.call_peaks_for_chromosome(chromosome, intervals)
            })
            .collect();
        
        #[cfg(not(feature = "parallel"))]
        let chromosome_peaks: Vec<Vec<Peak>> = self.coverage_data
            .iter()
            .map(|(chromosome, intervals)| {
                self.call_peaks_for_chromosome(chromosome, intervals)
            })
            .collect();

        let mut all_peaks: Vec<Peak> = chromosome_peaks.into_iter().flatten().collect();

        // Merge nearby peaks
        self.merge_nearby_peaks(all_peaks)
    }

    /// Call peaks for a single chromosome
    fn call_peaks_for_chromosome(&self, chromosome: &str, intervals: &[CoverageInterval]) -> Vec<Peak> {
        let mut peaks = Vec::new();
        let mut current_peak: Option<Peak> = None;

        // Process each position in the chromosome
        for interval in intervals {
            let signal = interval.coverage;
            let background = self.calculate_local_background(chromosome, interval.start);
            
            let fold_enrichment = if background > 0.0 {
                signal / background
            } else if signal > 0.0 {
                f64::INFINITY
            } else {
                0.0
            };

            // Check if this position meets peak criteria
            if fold_enrichment >= self.config.min_fold_enrichment {
                match current_peak.as_mut() {
                    Some(peak) => {
                        // Extend current peak
                        peak.end = interval.end;
                        if signal > peak.signal {
                            peak.signal = signal;
                            peak.background = background;
                            peak.fold_enrichment = fold_enrichment;
                        }
                    }
                    None => {
                        // Start new peak
                        current_peak = Some(Peak {
                            chromosome: chromosome.to_string(),
                            start: interval.start,
                            end: interval.end,
                            signal,
                            background,
                            fold_enrichment,
                        });
                    }
                }
            } else if let Some(peak) = current_peak.take() {
                // End current peak if it meets minimum length requirement
                if peak.end - peak.start >= self.config.min_peak_length {
                    peaks.push(peak);
                }
            }
        }

        // Don't forget the last peak
        if let Some(peak) = current_peak {
            if peak.end - peak.start >= self.config.min_peak_length {
                peaks.push(peak);
            }
        }

        peaks
    }

    /// Merge nearby peaks within merge_distance
    fn merge_nearby_peaks(&self, mut peaks: Vec<Peak>) -> Vec<Peak> {
        if peaks.is_empty() {
            return peaks;
        }

        // Sort peaks by chromosome and position
        peaks.sort_by(|a, b| {
            a.chromosome.cmp(&b.chromosome)
                .then_with(|| a.start.cmp(&b.start))
        });

        let mut merged_peaks = Vec::new();
        let mut current_peak = peaks[0].clone();

        for peak in peaks.into_iter().skip(1) {
            // Check if peaks are on same chromosome and close enough to merge
            if peak.chromosome == current_peak.chromosome 
                && peak.start <= current_peak.end + self.config.merge_distance {
                // Merge peaks
                current_peak.end = peak.end.max(current_peak.end);
                if peak.signal > current_peak.signal {
                    current_peak.signal = peak.signal;
                    current_peak.background = peak.background;
                    current_peak.fold_enrichment = peak.fold_enrichment;
                }
            } else {
                // Save current peak and start new one
                merged_peaks.push(current_peak);
                current_peak = peak;
            }
        }

        merged_peaks.push(current_peak);
        merged_peaks
    }

    /// Write peaks to BED format file
    pub fn write_peaks_to_bed<P: AsRef<Path>>(&self, peaks: &[Peak], output_path: P) -> Result<(), Box<dyn std::error::Error>> {
        let mut file = File::create(output_path)?;
        
        // Write BED header
        writeln!(file, "track name=\"ChIP-Seq Peaks\" description=\"Peaks called with fold enrichment >= {:.1}\"", 
                self.config.min_fold_enrichment)?;

        for peak in peaks {
            // BED format: chromosome start end name score strand
            writeln!(
                file,
                "{}\t{}\t{}\tpeak\t{:.0}\t.",
                peak.chromosome,
                peak.start,
                peak.end,
                (peak.fold_enrichment * 100.0).min(1000.0) // Score capped at 1000
            )?;
        }

        Ok(())
    }

    /// Get summary statistics
    pub fn get_summary(&self, peaks: &[Peak]) -> PeakCallingSummary {
        if peaks.is_empty() {
            return PeakCallingSummary::default();
        }

        let total_peaks = peaks.len();
        let total_coverage: u32 = peaks.iter().map(|p| p.end - p.start).sum();
        let avg_length = total_coverage as f64 / total_peaks as f64;
        let max_enrichment = peaks.iter().map(|p| p.fold_enrichment).fold(0.0, f64::max);
        let avg_enrichment = peaks.iter().map(|p| p.fold_enrichment).sum::<f64>() / total_peaks as f64;

        // Count peaks per chromosome
        let mut peaks_per_chromosome = BTreeMap::new();
        for peak in peaks {
            *peaks_per_chromosome.entry(peak.chromosome.clone()).or_insert(0) += 1;
        }

        PeakCallingSummary {
            total_peaks,
            total_coverage,
            average_peak_length: avg_length,
            max_fold_enrichment: max_enrichment,
            average_fold_enrichment: avg_enrichment,
            peaks_per_chromosome,
        }
    }
}

/// Summary statistics for peak calling results
#[derive(Debug, Default)]
pub struct PeakCallingSummary {
    pub total_peaks: usize,
    pub total_coverage: u32,
    pub average_peak_length: f64,
    pub max_fold_enrichment: f64,
    pub average_fold_enrichment: f64,
    pub peaks_per_chromosome: BTreeMap<String, usize>,
}

impl std::fmt::Display for PeakCallingSummary {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Peak Calling Summary:")?;
        writeln!(f, "  Total peaks: {}", self.total_peaks)?;
        writeln!(f, "  Total coverage: {} bp", self.total_coverage)?;
        writeln!(f, "  Average peak length: {:.1} bp", self.average_peak_length)?;
        writeln!(f, "  Maximum fold enrichment: {:.2}", self.max_fold_enrichment)?;
        writeln!(f, "  Average fold enrichment: {:.2}", self.average_fold_enrichment)?;
        writeln!(f, "  Peaks per chromosome:")?;
        for (chr, count) in &self.peaks_per_chromosome {
            writeln!(f, "    {}: {}", chr, count)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_peak_calling_basic() {
        let mut caller = PeakCaller::new();
        
        // Add test coverage data
        caller.add_coverage_interval(CoverageInterval {
            chromosome: "chr1".to_string(),
            start: 1000,
            end: 1100,
            coverage: 10.0,
        });
        
        caller.add_coverage_interval(CoverageInterval {
            chromosome: "chr1".to_string(),
            start: 1100,
            end: 1200,
            coverage: 50.0, // This should be a peak
        });
        
        caller.add_coverage_interval(CoverageInterval {
            chromosome: "chr1".to_string(),
            start: 1200,
            end: 1300,
            coverage: 8.0,
        });

        let peaks = caller.call_peaks();
        assert!(!peaks.is_empty(), "Should find at least one peak");
        
        let peak = &peaks[0];
        assert_eq!(peak.chromosome, "chr1");
        assert!(peak.fold_enrichment >= 2.0);
    }

    #[test]
    fn test_bedgraph_loading() -> Result<(), Box<dyn std::error::Error>> {
        let mut temp_file = NamedTempFile::new()?;
        writeln!(temp_file, "chr1\t1000\t1100\t5.5")?;
        writeln!(temp_file, "chr1\t1100\t1200\t25.0")?;
        writeln!(temp_file, "chr2\t2000\t2100\t12.0")?;
        
        let mut caller = PeakCaller::new();
        caller.load_coverage_from_bedgraph(temp_file.path())?;
        
        assert_eq!(caller.coverage_data.len(), 2);
        assert!(caller.coverage_data.contains_key("chr1"));
        assert!(caller.coverage_data.contains_key("chr2"));
        
        Ok(())
    }

    #[test]
    fn test_local_background_calculation() {
        let mut caller = PeakCaller::with_config(PeakCallingConfig {
            background_window: 1000,
            ..Default::default()
        });
        
        // Add uniform background coverage
        for i in 0..10 {
            caller.add_coverage_interval(CoverageInterval {
                chromosome: "chr1".to_string(),
                start: i * 100,
                end: (i + 1) * 100,
                coverage: 5.0,
            });
        }
        
        let background = caller.calculate_local_background("chr1", 500);
        assert!((background - 5.0).abs() < 0.1, "Background should be ~5.0, got {}", background);
    }

    #[test] 
    fn test_peak_merging() {
        let mut caller = PeakCaller::with_config(PeakCallingConfig {
            merge_distance: 150,
            min_peak_length: 50,
            ..Default::default()
        });
        
        // Create two close peaks that should be merged
        caller.add_coverage_interval(CoverageInterval {
            chromosome: "chr1".to_string(),
            start: 1000,
            end: 1100,
            coverage: 20.0,
        });
        
        caller.add_coverage_interval(CoverageInterval {
            chromosome: "chr1".to_string(),
            start: 1200, // 100bp gap, within merge distance
            end: 1300,
            coverage: 25.0,
        });
        
        let peaks = caller.call_peaks();
        // Depending on background, these might merge into 1 peak
        assert!(!peaks.is_empty());
    }

    #[test]
    fn test_bed_output() -> Result<(), Box<dyn std::error::Error>> {
        let caller = PeakCaller::new();
        let peaks = vec![
            Peak {
                chromosome: "chr1".to_string(),
                start: 1000,
                end: 1200,
                signal: 50.0,
                background: 10.0,
                fold_enrichment: 5.0,
            }
        ];
        
        let temp_file = NamedTempFile::new()?;
        caller.write_peaks_to_bed(&peaks, temp_file.path())?;
        
        let content = std::fs::read_to_string(temp_file.path())?;
        assert!(content.contains("chr1\t1000\t1200"));
        assert!(content.contains("track name"));
        
        Ok(())
    }
}

// Example usage
fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Create peak caller with custom configuration
    let config = PeakCallingConfig {
        min_fold_enrichment: 3.0,
        background_window: 5000,
        min_peak_length: 200,
        merge_distance: 300,
    };
    
    let mut caller = PeakCaller::with_config(config);
    
    // Load coverage data from BedGraph file
    // caller.load_coverage_from_bedgraph("input_coverage.bedgraph")?;
    
    // Or add data programmatically for demonstration
    caller.add_coverage_interval(CoverageInterval {
        chromosome: "chr1".to_string(),
        start: 1000,
        end: 1200,
        coverage: 5.0,
    });
    
    caller.add_coverage_interval(CoverageInterval {
        chromosome: "chr1".to_string(),
        start: 1200,
        end: 1400,
        coverage: 30.0, // Strong signal
    });
    
    caller.add_coverage_interval(CoverageInterval {
        chromosome: "chr1".to_string(),
        start: 1400,
        end: 1600,
        coverage: 25.0, // Continued signal
    });
    
    caller.add_coverage_interval(CoverageInterval {
        chromosome: "chr1".to_string(),
        start: 1600,
        end: 1800,
        coverage: 4.0,
    });
    
    // Call peaks
    let peaks = caller.call_peaks();
    
    // Print summary
    let summary = caller.get_summary(&peaks);
    println!("{}", summary);
    
    // Write peaks to BED file
    caller.write_peaks_to_bed(&peaks, "output_peaks.bed")?;
    
    println!("Peak calling completed! Found {} peaks.", peaks.len());
    for (i, peak) in peaks.iter().enumerate() {
        println!("Peak {}: {}:{}-{} (fold enrichment: {:.2})", 
                i + 1, peak.chromosome, peak.start, peak.end, peak.fold_enrichment);
    }
    
    Ok(())
}