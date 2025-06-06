use std::collections::HashMap;
use std::path::PathBuf;
use std::fs::File;
use std::io::Write;
use anyhow::{Context, Result};
use clap::Parser;
use log::info;
use rayon::prelude::*;
use rust_htslib::bam::{self, Read};
use serde::{Deserialize, Serialize};

#[derive(Parser, Debug)]
#[command(author, version, about = "Simple peak caller for QC pipeline")]
struct Args {
    /// Input BAM file
    input: PathBuf,
    
    /// Output BED file (default: stdout)
    #[arg(short, long)]
    output: Option<PathBuf>,
    
    /// Window size for peak calling (bp)
    #[arg(short, long, default_value_t = 1000)]
    window_size: u64,
    
    /// Minimum coverage threshold
    #[arg(short, long, default_value_t = 5.0)]
    min_coverage: f64,
    
    /// Number of threads
    #[arg(short, long, default_value_t = 4)]
    threads: usize,
}

#[derive(Debug, Serialize, Deserialize)]
struct SimplePeak {
    chrom: String,
    start: u64,
    end: u64,
    coverage: f64,
}

fn main() -> Result<()> {
    let args = Args::parse();
    
    env_logger::init();
    
    info!("Starting simple peak calling...");
    info!("Input: {:?}", args.input);
    info!("Window size: {} bp", args.window_size);
    info!("Min coverage: {}", args.min_coverage);
    
    // Set thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .context("Failed to initialize thread pool")?;
    
    // Read BAM file and calculate coverage
    let mut bam_reader = bam::Reader::from_path(&args.input)
        .with_context(|| format!("Failed to open BAM file: {:?}", args.input))?;
    
    let header = bam_reader.header().clone();
    
    // Collect read positions
    let mut coverage_data: HashMap<String, HashMap<u64, u64>> = HashMap::new();
    
    info!("Processing reads...");
    for result in bam_reader.records() {
        let record = result?;
        
        if record.is_unmapped() || record.is_secondary() || record.is_duplicate() {
            continue;
        }
        
        let tid = record.tid();
        if tid < 0 {
            continue;
        }
        
        let chrom_name = String::from_utf8_lossy(header.tid2name(tid as u32)?).to_string();
        let start = record.pos() as u64;
        let window = (start / args.window_size) * args.window_size;
        
        let chrom_map = coverage_data.entry(chrom_name).or_insert_with(HashMap::new);
        *chrom_map.entry(window).or_insert(0) += 1;
    }
    
    info!("Calling peaks...");
    
    // Find peaks
    let mut peaks = Vec::new();
    for (chrom, windows) in coverage_data {
        for (window_start, count) in windows {
            let coverage = count as f64;
            if coverage >= args.min_coverage {
                peaks.push(SimplePeak {
                    chrom: chrom.clone(),
                    start: window_start,
                    end: window_start + args.window_size,
                    coverage,
                });
            }
        }
    }
    
    // Sort peaks by coverage (descending)
    peaks.sort_by(|a, b| b.coverage.partial_cmp(&a.coverage).unwrap());
    
    info!("Found {} peaks", peaks.len());
    
    // Write output
    let output_content = peaks
        .iter()
        .map(|peak| format!("{}\t{}\t{}\tpeak\t{}", 
                           peak.chrom, peak.start, peak.end, peak.coverage as u32))
        .collect::<Vec<_>>()
        .join("\n");
    
    if let Some(output_path) = args.output {
        let mut file = File::create(output_path)?;
        writeln!(file, "{}", output_content)?;
    } else {
        println!("{}", output_content);
    }
    
    info!("Peak calling completed");
    Ok(())
}