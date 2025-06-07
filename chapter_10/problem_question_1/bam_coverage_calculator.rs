use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use rayon::prelude::*;
use rust_htslib::{bam, bam::Read};
use clap::{Arg, Command};

#[derive(Debug, Clone)]
struct CoverageInterval {
    start: u64,
    end: u64,
    depth: u32,
}

#[derive(Debug)]
struct ChromosomeCoverage {
    chromosome: String,
    intervals: Vec<CoverageInterval>,
}

impl ChromosomeCoverage {
    fn new(chromosome: String) -> Self {
        Self {
            chromosome,
            intervals: Vec::new(),
        }
    }

    fn add_interval(&mut self, start: u64, end: u64, depth: u32) {
        self.intervals.push(CoverageInterval { start, end, depth });
    }

    fn merge_intervals(&mut self) {
        if self.intervals.is_empty() {
            return;
        }

        // Sort intervals by start position
        self.intervals.sort_by_key(|interval| interval.start);

        let mut merged = Vec::new();
        let mut current = self.intervals[0].clone();

        for interval in self.intervals.iter().skip(1) {
            if interval.start <= current.end && interval.depth == current.depth {
                // Merge adjacent intervals with same depth
                current.end = current.end.max(interval.end);
            } else {
                merged.push(current);
                current = interval.clone();
            }
        }
        merged.push(current);
        self.intervals = merged;
    }
}

struct CoverageCalculator {
    input_path: String,
    output_path: String,
    min_depth: u32,
}

impl CoverageCalculator {
    fn new(input_path: String, output_path: String, min_depth: u32) -> Self {
        Self {
            input_path,
            output_path,
            min_depth,
        }
    }

    fn calculate_coverage(&self) -> Result<(), Box<dyn std::error::Error>> {
        println!("Reading BAM file: {}", self.input_path);
        
        // Open BAM file
        let mut bam_reader = bam::Reader::from_path(&self.input_path)?;
        let header = bam_reader.header().clone();
        
        // Get chromosome names from header
        let chromosomes: Vec<String> = (0..header.target_count())
            .map(|i| String::from_utf8_lossy(header.target_names()[i as usize]).to_string())
            .collect();

        println!("Found {} chromosomes", chromosomes.len());

        // Process each chromosome separately to manage memory efficiently
        let mut all_coverage: Vec<ChromosomeCoverage> = Vec::new();

        for (chr_idx, chr_name) in chromosomes.iter().enumerate() {
            println!("Processing chromosome: {}", chr_name);
            
            let chr_coverage = self.process_chromosome(&mut bam_reader, chr_idx as u32, chr_name)?;
            if !chr_coverage.intervals.is_empty() {
                all_coverage.push(chr_coverage);
            }
        }

        // Sort chromosomes and write output
        self.write_coverage_bed(&all_coverage)?;
        
        Ok(())
    }

    fn process_chromosome(
        &self,
        bam_reader: &mut bam::Reader,
        chr_idx: u32,
        chr_name: &str,
    ) -> Result<ChromosomeCoverage, Box<dyn std::error::Error>> {
        
        // Use a more memory-efficient approach with a sliding window
        let mut coverage_map: HashMap<u64, u32> = HashMap::new();
        let mut record = bam::Record::new();

        // Fetch reads for this chromosome
        bam_reader.fetch(bam::FetchDefinition::RegionIndex(chr_idx))?;

        while let Some(result) = bam_reader.records().next() {
            let record = result?;
            
            // Skip unmapped reads
            if record.is_unmapped() {
                continue;
            }

            let start = record.pos() as u64;
            let end = record.cigar_cached().unwrap().end_pos() as u64;

            // Increment coverage for each position covered by this read
            for pos in start..end {
                *coverage_map.entry(pos).or_insert(0) += 1;
            }
        }

        // Convert coverage map to intervals
        let mut chr_coverage = ChromosomeCoverage::new(chr_name.to_string());
        
        if coverage_map.is_empty() {
            return Ok(chr_coverage);
        }

        // Sort positions and create intervals
        let mut positions: Vec<u64> = coverage_map.keys().cloned().collect();
        positions.sort();

        let mut current_start = positions[0];
        let mut current_end = positions[0] + 1;
        let mut current_depth = coverage_map[&positions[0]];

        for &pos in positions.iter().skip(1) {
            let depth = coverage_map[&pos];
            
            if pos == current_end && depth == current_depth {
                // Extend current interval
                current_end = pos + 1;
            } else {
                // End current interval and start new one
                if current_depth >= self.min_depth {
                    chr_coverage.add_interval(current_start, current_end, current_depth);
                }
                current_start = pos;
                current_end = pos + 1;
                current_depth = depth;
            }
        }

        // Add the last interval
        if current_depth >= self.min_depth {
            chr_coverage.add_interval(current_start, current_end, current_depth);
        }

        chr_coverage.merge_intervals();
        Ok(chr_coverage)
    }

    fn write_coverage_bed(&self, coverage_data: &[ChromosomeCoverage]) -> Result<(), Box<dyn std::error::Error>> {
        println!("Writing coverage to: {}", self.output_path);
        
        let file = File::create(&self.output_path)?;
        let mut writer = BufWriter::new(file);

        // Write BED header
        writeln!(writer, "track type=bedGraph name=\"Coverage Track\" description=\"BAM Coverage\"")?;

        for chr_coverage in coverage_data {
            for interval in &chr_coverage.intervals {
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}",
                    chr_coverage.chromosome,
                    interval.start,
                    interval.end,
                    interval.depth
                )?;
            }
        }

        writer.flush()?;
        println!("Coverage calculation complete!");
        Ok(())
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let matches = Command::new("BAM Coverage Calculator")
        .version("1.0")
        .author("Bioinformatics Tool")
        .about("Calculates coverage track from BAM files")
        .arg(
            Arg::new("input")
                .short('i')
                .long("input")
                .value_name("BAM_FILE")
                .help("Input BAM file path")
                .required(true)
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_name("BED_FILE")
                .help("Output BED file path")
                .required(true)
        )
        .arg(
            Arg::new("min-depth")
                .short('d')
                .long("min-depth")
                .value_name("DEPTH")
                .help("Minimum coverage depth to report")
                .default_value("1")
        )
        .get_matches();

    let input_path = matches.get_one::<String>("input").unwrap().to_string();
    let output_path = matches.get_one::<String>("output").unwrap().to_string();
    let min_depth: u32 = matches.get_one::<String>("min-depth").unwrap().parse()?;

    // Validate input file exists
    if !Path::new(&input_path).exists() {
        eprintln!("Error: Input BAM file '{}' does not exist", input_path);
        std::process::exit(1);
    }

    let calculator = CoverageCalculator::new(input_path, output_path, min_depth);
    calculator.calculate_coverage()?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_interval_merging() {
        let mut chr_coverage = ChromosomeCoverage::new("chr1".to_string());
        
        // Add overlapping intervals with same depth
        chr_coverage.add_interval(100, 200, 5);
        chr_coverage.add_interval(150, 250, 5);
        chr_coverage.add_interval(300, 400, 3);
        
        chr_coverage.merge_intervals();
        
        assert_eq!(chr_coverage.intervals.len(), 2);
        assert_eq!(chr_coverage.intervals[0].start, 100);
        assert_eq!(chr_coverage.intervals[0].end, 250);
        assert_eq!(chr_coverage.intervals[0].depth, 5);
    }

    #[test]
    fn test_coverage_interval_creation() {
        let interval = CoverageInterval {
            start: 1000,
            end: 2000,
            depth: 10,
        };
        
        assert_eq!(interval.start, 1000);
        assert_eq!(interval.end, 2000);
        assert_eq!(interval.depth, 10);
    }
}