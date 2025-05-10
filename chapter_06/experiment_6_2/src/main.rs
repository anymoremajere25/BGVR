use anyhow::{Context, Result};
use clap::Parser;
use env_logger;
use log::{error, info};
use rayon::prelude::*;
use rust_htslib::bam::{self, Read};
use std::path::Path;
use std::error::Error;

/// A command-line interface for the coverage calculator.
#[derive(Parser, Debug)]
#[command(name = "coverage_tool")]
#[command(author, version, about = "A robust coverage calculator using rust-htslib and rayon")]
struct Cli {
    /// Path to the input BAM file (must be indexed)
    #[arg(long, value_name = "BAM_FILE")]
    bam: String,

    /// Genomic regions to calculate coverage for, e.g., 1:1-50000
    #[arg(long, value_name = "REGIONS", required = true, num_args = 1..)]
    region: Vec<String>,
}

/// Calculate coverage for multiple regions in parallel using an indexed BAM file.
/// Returns a vector of (region, coverage_sum) tuples for each region.
pub fn parallel_coverage(bam_path: &str, regions: &[String]) -> Result<Vec<(String, u64)>> {
    let results: Vec<(String, u64)> = regions
        .par_iter()
        .map(|region| {
            // Create a new IndexedReader for each region to ensure thread safety
            let mut reader = bam::IndexedReader::from_path(&Path::new(bam_path))
                .with_context(|| format!("Failed to open indexed BAM file: {bam_path}"))?;

            // Fetch the specified region from the BAM/CRAM file
            reader
                .fetch(region)
                .with_context(|| format!("Failed to fetch region: {region}"))?;

            // Accumulate coverage by summing the lengths of all read sequences in the region
            let mut coverage_sum = 0u64;
            for rec in reader.records() {
                let record = rec.with_context(|| {
                    format!("Error reading record in region {region} from file {bam_path}")
                })?;
                coverage_sum += record.seq().len() as u64;
            }
            Ok((region.clone(), coverage_sum))
        })
        .collect::<Result<Vec<(String, u64)>>>()?;

    Ok(results)
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();

    // Log the file and regions that will be analyzed
    info!("Starting coverage analysis on file: {}", cli.bam);
    info!("Regions to analyze: {:?}", cli.region);

    // Calculate coverage in parallel using our function
    let coverage_info = match parallel_coverage(&cli.bam, &cli.region) {
        Ok(info) => info,
        Err(e) => {
            error!("Coverage calculation failed: {:?}", e);
            return Err(e);
        }
    };

    // Print or otherwise process the coverage results
    for (region, sum) in coverage_info {
        println!("Region: {}, Coverage: {}", region, sum);
    }

    info!("Coverage analysis completed successfully.");
    Ok(())
}
