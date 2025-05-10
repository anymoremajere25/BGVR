use anyhow::{Context, Result};
use clap::Parser;
use rust_htslib::bcf::{Format, Read, Reader, Writer, Header};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::sync::{Arc, Mutex};
use std::time::Instant;
use rayon::prelude::*;

/// Command line arguments
#[derive(Parser, Debug)]
#[command(author, version, about = "Merge VCF files into a multi-sample VCF")]
struct Args {
    #[arg(short, long)]
    vcfs: Vec<String>,

    #[arg(long)]
    vcf_list: Option<String>, // Optional file with list of VCF paths

    #[arg(short, long)]
    out: String,

    #[arg(short, long, default_value = "1")]
    threads: usize,

    #[arg(short, long)]
    skip_malformed: bool,

    #[arg(long)]
    dry_run: bool,

    #[arg(long, default_value = "bcf")]
    format: String,
}

fn load_vcf_files(args: &Args) -> Result<Vec<String>> {
    let mut vcf_files = args.vcfs.clone();

    if let Some(list_path) = &args.vcf_list {
        let file = File::open(list_path).with_context(|| "Failed to open VCF list file")?;
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let path = line?.trim().to_string();
            if !path.is_empty() {
                vcf_files.push(path);
            }
        }
    }

    Ok(vcf_files)
}

fn process_vcf_file(path: &str, writer: &Arc<Mutex<Writer>>) -> Result<()> {
    let mut reader = Reader::from_path(path)
        .with_context(|| format!("Failed to open VCF file: {}", path))?;

    for rec in reader.records() {
        match rec {
            Ok(rec) => writer.lock().unwrap().write(&rec)?,
            Err(err) => println!("Failed to read record in {}: {}", path, err),
        }
    }

    Ok(())
}

fn merge_vcfs(vcf_files: Vec<String>, output_path: &str, threads: usize, format: &str, dry_run: bool) -> Result<()> {
    let start_time = Instant::now();

    if dry_run {
        println!("Dry run complete — no output written.");
        return Ok(());
    }

    let mut reader = Reader::from_path(&vcf_files[0])?;
    let header = Header::from_template(reader.header());

    // Match string to rust_htslib format enum
    let format_enum = match format.to_lowercase().as_str() {
        "vcf" => Format::Vcf,
        "bcf" => Format::Bcf,
        _ => return Err(anyhow::anyhow!("Unsupported format: {}", format)),
    };

    let writer = Writer::from_path(output_path, &header, false, format_enum)?;
    let writer = Arc::new(Mutex::new(writer));

    if threads > 1 {
        vcf_files.par_iter().try_for_each(|path| process_vcf_file(path, &writer))?;
    } else {
        for path in &vcf_files {
            process_vcf_file(path, &writer)?;
        }
    }

    println!("Merge completed in {} ms", start_time.elapsed().as_millis());
    Ok(())
}

fn main() -> Result<()> {
    let args = Args::parse();
    let vcf_files = load_vcf_files(&args)?;
    merge_vcfs(vcf_files, &args.out, args.threads, &args.format, args.dry_run)
}