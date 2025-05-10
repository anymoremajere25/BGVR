use anyhow::{Context, Result};
use clap::Parser;
use memmap2::{Mmap, MmapOptions};
use rayon::prelude::*;
use serde::Serialize;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::str;
use std::time::Instant;

#[derive(Parser, Debug)]
#[clap(about = "High-performance genomic analysis tool")]
struct Args {
    #[clap(short = 'f', long, default_value = "reference.fasta")]
    reference: PathBuf,

    #[clap(short = 'r', long)]
    region: Option<String>,

    #[clap(short, long)]
    output: Option<PathBuf>,

    #[clap(short, long)]
    threads: Option<usize>,

    #[clap(short, long)]
    verbose: bool,
}

#[derive(Debug, Serialize)]
struct AnalysisResult {
    region: Option<String>,
    gc_content: f64,
    sequence_length: usize,
}

fn main() -> Result<()> {
    let args = Args::parse();

    env_logger::Builder::from_default_env()
        .filter_level(if args.verbose { log::LevelFilter::Debug } else { log::LevelFilter::Info })
        .init();

    if let Some(threads) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()?;
    }

    let start = Instant::now();
    println!("Processing file: {:?}", args.reference);

    let file = File::open(&args.reference)?;
    let mmap = unsafe { MmapOptions::new().map(&file)? };
    let result = process_fasta(&mmap, args.region.as_deref())?;

    match args.output {
        Some(path) => write_results_to_file(&result, &path)?,
        None => println!("{}", serde_json::to_string_pretty(&result)?),
    }

    println!("Analysis completed in {:?}", start.elapsed());
    Ok(())
}

fn process_fasta(mmap: &Mmap, region_str: Option<&str>) -> Result<AnalysisResult> {
    let contents = str::from_utf8(mmap).context("Invalid UTF-8 in FASTA")?;
    let lines: Vec<&str> = contents.lines().collect();

    // Split header and sequence lines
    let (headers, sequences): (Vec<&str>, Vec<&str>) = lines.par_iter()
        .partition(|line| line.starts_with('>'));

    println!("Found {} sequences in FASTA", headers.len());

    // Calculate GC content in parallel
    let (gc_count, total_bases) = sequences.par_iter()
        .map(|seq| {
            let seq_bytes = seq.as_bytes();
            let gc = seq_bytes.iter().filter(|&&b| b == b'G' || b == b'g' || b == b'C' || b == b'c').count();
            (gc, seq_bytes.len())
        })
        .reduce(|| (0, 0), |(gc1, len1), (gc2, len2)| (gc1 + gc2, len1 + len2));

    let gc_content = if total_bases > 0 {
        gc_count as f64 / total_bases as f64
    } else {
        0.0
    };

    Ok(AnalysisResult {
        region: region_str.map(String::from),
        gc_content,
        sequence_length: total_bases,
    })
}

fn write_results_to_file(results: &AnalysisResult, path: &Path) -> Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    writer.write_all(serde_json::to_string_pretty(results)?.as_bytes())?;
    println!("Results written to {:?}", path);
    Ok(())
}