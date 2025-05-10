use std::path::Path;
use anyhow::{Context, Result};
use clap::Parser;
use env_logger;
use log::{info, warn, error};
use rayon::prelude::*;
use rust_htslib::bcf::{Format, Header, Reader, Record, Writer};
use rust_htslib::bcf::Read; // Needed for header() and records()

/// Command-line arguments for the filtering tool
#[derive(Parser, Debug)]
#[command(name = "bcf_filter_tool", version, about = "A robust BCF filtering tool using rust-htslib and rayon")]
struct Cli {
    /// Input BCF/VCF file (must be indexed)
    #[arg(long, value_name = "INPUT_BCF", short = 'i')]
    input: String,

    /// Output file path where filtered records should be written
    #[arg(long, value_name = "OUTPUT_BCF", short = 'o')]
    output: String,

    /// Minimum QUAL score to retain
    #[arg(long, value_name = "MIN_QUAL", default_value = "30.0")]
    min_qual: f64,

    /// Minimum average depth (DP) to retain
    #[arg(long, value_name = "MIN_DEPTH", default_value = "10")]
    min_depth: i32,

    /// Number of records to process at a time
    #[arg(long, value_name = "CHUNK_SIZE", default_value = "10000")]
    chunk_size: usize,
}

/// Parameters for filtering variant records
#[derive(Debug)]
struct FilterParams {
    min_qual: f64,
    min_depth: i32,
}

/// Filters a batch of records in parallel based on a minimum QUAL and DP threshold
fn filter_records(records: Vec<Record>, params: &FilterParams) -> Vec<Record> {
    records
        .into_par_iter()
        .filter(|record| {
            // Filter on QUAL
            if record.qual() < params.min_qual as f32 {
                return false;
            }

            // Filter on average DP
            if let Ok(depths) = record.format(b"DP").integer() {
                let dps: Vec<i32> = depths.iter().filter_map(|arr| arr.get(0).copied()).collect();
                if !dps.is_empty() {
                    let avg_dp = dps.iter().map(|&d| d as f64).sum::<f64>() / (dps.len() as f64);
                    if avg_dp < params.min_depth as f64 {
                        return false;
                    }
                } else {
                    warn!("Record found without DP field. Retaining by default.");
                }
            }
            true
        })
        .collect()
}

/// Filters a BCF file in parallel and writes the result to another file.
fn parallel_filter_vcf(
    input_path: &str,
    output_path: &str,
    params: &FilterParams,
    chunk_size: usize,
) -> Result<()> {
    let mut reader = Reader::from_path(Path::new(input_path))
        .with_context(|| format!("Failed to open BCF file: {}", input_path))?;

    // Convert header from HeaderView to Header
    let header = Header::from_template(reader.header());
    let mut writer = Writer::from_path(output_path, &header, false, Format::Bcf)
        .with_context(|| format!("Failed to create writer for output: {}", output_path))?;

    let mut record_batch = Vec::with_capacity(chunk_size);
    let mut total_processed = 0usize;
    let mut total_retained = 0usize;

    for result in reader.records() {
        let record = result.with_context(|| {
            format!("Error reading a record from file: {}", input_path)
        })?;

        record_batch.push(record);
        if record_batch.len() == chunk_size {
            let filtered_records = filter_records(record_batch.drain(..).collect(), params);
            total_retained += filtered_records.len();
            for rec in filtered_records {
                writer.write(&rec).with_context(|| "Error writing filtered record")?;
            }
            total_processed += chunk_size;
        }
    }

    if !record_batch.is_empty() {
        let filtered_records = filter_records(record_batch, params);
        total_retained += filtered_records.len();
        for rec in &filtered_records {
            writer.write(rec).with_context(|| "Error writing final batch of records")?;
        }
        total_processed += filtered_records.len();
    }
    
    info!(
        "Filtering completed. Processed {} records in total. Retained {} records after filtering.",
        total_processed, total_retained
    );

    Ok(())
}

fn main() -> Result<()> {
    env_logger::init();

    let cli = Cli::parse();
    let params = FilterParams {
        min_qual: cli.min_qual,
        min_depth: cli.min_depth,
    };

    info!("Starting parallel BCF filtering with parameters: {:?}", params);

    match parallel_filter_vcf(&cli.input, &cli.output, &params, cli.chunk_size) {
        Ok(_) => info!("Filtering completed successfully."),
        Err(err) => {
            error!("Filtering failed with error: {:?}", err);
            return Err(err);
        }
    }

    Ok(())
}