// src/bin/bam_indexer.rs
use clap::Parser;
use std::process::Command;

#[derive(Parser)]
#[command(name = "bam_indexer")]
#[command(about = "Index BAM files for efficient random access")]
pub struct Args {
    #[arg(long, help = "Input BAM file")]
    input: String,
    
    #[arg(long, help = "Output BAI index file")]
    output: String,
    
    #[arg(long, default_value = "4", help = "Number of threads")]
    threads: usize,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    
    println!("Indexing BAM file: {}", args.input);
    println!("Output index: {}", args.output);
    
    // Use samtools to index the BAM file
    let mut cmd = Command::new("samtools");
    cmd.args(&[
        "index",
        "-@", &args.threads.to_string(),
        &args.input,
        &args.output,
    ]);
    
    let output = cmd.output()?;
    
    if output.status.success() {
        println!("BAM file indexed successfully");
        Ok(())
    } else {
        let error_msg = String::from_utf8_lossy(&output.stderr);
        eprintln!("BAM indexing failed: {}", error_msg);
        Err(format!("BAM indexing failed with exit code: {}", output.status).into())
    }
}