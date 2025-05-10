use anyhow::{Context, Result};
use clap::Parser;
use rayon::prelude::*;
use rust_htslib::bam::{IndexedReader, Read};
use std::path::Path;

/// A structure that keeps track of total and mapped read counts.
#[derive(Debug)]
struct ReadCounts {
    total_reads: u64,
    mapped_reads: u64,
}

/// Command-line interface to set up parameters for parallel read counting.
#[derive(Parser, Debug)]
#[command(name = "bam_read_counter")]
#[command(author, version, about = "A robust tool to count total and mapped reads across multiple regions")]
struct Cli {
    /// Path to the indexed BAM file (must have a corresponding BAI file)
    #[arg(long, value_name = "BAM_FILE")]
    bam: String,

    /// One or more regions (e.g., 'chr1:1-1000000') for counting reads
    #[arg(long, value_name = "REGION", required = true, num_args = 1..)]
    region: Vec<String>,
}

/// Processes a single region by fetching reads from the BAM file and counting total vs. mapped reads.
fn process_bam_chunk(bam_path: &str, region: &str) -> Result<ReadCounts> {
    // Open the BAM file with an indexed reader
    let mut reader = IndexedReader::from_path(Path::new(bam_path))
        .with_context(|| format!("Failed to open indexed BAM file: {}", bam_path))?;
    
    // Restrict reading to the specified region
    reader.fetch(region)
        .with_context(|| format!("Failed to fetch region {region}"))?;

    let mut counts = ReadCounts { total_reads: 0, mapped_reads: 0 };
    // Iterate over each record in the region
    for result in reader.records() {
        let record = result.with_context(|| format!("Error reading record in region {region}"))?;
        counts.total_reads += 1;
        if !record.is_unmapped() {
            counts.mapped_reads += 1;
        }
    }
    Ok(counts)
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    println!("Starting parallel read counting on file: {}", cli.bam);
    println!("Processing regions: {:?}", cli.region);

    // Parallel execution over the specified regions
    let results: Vec<Result<ReadCounts>> = cli.region
        .par_iter()
        .map(|region| process_bam_chunk(&cli.bam, region))
        .collect();

    // Aggregate results and handle any errors encountered in parallel tasks
    let mut total_reads: u64 = 0;
    let mut mapped_reads: u64 = 0;
    let mut successful_regions = 0;
    for (region, res) in cli.region.iter().zip(results.into_iter()) {
        match res {
            Ok(counts) => {
                total_reads += counts.total_reads;
                mapped_reads += counts.mapped_reads;
                successful_regions += 1;
            }
            Err(e) => {
                eprintln!("Failed to process region {}: {:?}", region, e);
            }
        }
    }

    println!("Successfully processed {} out of {} regions.", successful_regions, cli.region.len());
    println!("Total reads: {}, Mapped reads: {}", total_reads, mapped_reads);

    Ok(())
}
