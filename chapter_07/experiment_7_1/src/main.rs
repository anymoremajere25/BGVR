use noodles::bam;
use std::fs::File;
use std::path::Path;
use std::error::Error;
use clap::{Arg, Command};
use rayon::prelude::*;
use log::error;
use anyhow::Result;

// Function to read BAM file and process it
fn read_bam_file(file_path: &str) -> Result<(), Box<dyn Error>> {
    println!("Processing BAM file: {}", file_path);
    let path = Path::new(file_path);
    
    // Open the file manually
    let file = File::open(path)?;
    
    // Create the reader using noodles::bam
    let mut reader = bam::Reader::new(file);
    
    // Read the header first
    let header = reader.read_header()?;
    println!("Successfully read BAM header");
    
    // Read and process the reference sequences
    let reference_sequences = reader.read_reference_sequences()?;
    println!("Reference sequences count: {}", reference_sequences.len());
    
    // Print some reference sequence info
    for (name, seq) in reference_sequences.iter().take(2) {
        println!("  Reference: {} (length: {})", name, seq.len());
    }
    
    // Use a more careful approach to read records
    let mut record_count = 0;
    
    // Create a new record outside the loop to reuse
    let mut record = bam::Record::default();
    
    // Read records
    while reader.read_record(&mut record)? > 0 {
        record_count += 1;
        // Process specific fields from the record
        if record_count <= 5 {  // Only print first 5 records to avoid overwhelming output
            // Access basic record information
            let position = record.position();
            let mapping_quality = record.mapping_quality().unwrap_or(0);
            let cigar = record.cigar();
            
            println!("Record #{}: Position={:?}, MAPQ={}, CIGAR={:?}", 
                     record_count, position, mapping_quality, cigar);
        }
    }
    
    println!("Total records processed from {}: {}", file_path, record_count);
    Ok(())
}

// Function to parse VCF file
fn read_vcf_file(vcf_file: &str) -> Result<(), Box<dyn Error>> {
    // Since VCF feature isn't included, we'll just acknowledge the file
    println!("VCF file noted: {} (VCF processing not enabled)", vcf_file);
    println!("Add 'vcf' to noodles features in Cargo.toml to enable VCF processing");
    
    Ok(())
}

// Initialize logging with specified level
fn init_logging() -> Result<(), Box<dyn Error>> {
    env_logger::Builder::new()
        .filter_level(log::LevelFilter::Info)
        .init();
    
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    // Initialize logging
    init_logging()?;
    
    println!("Starting BAM processing application...");
    
    // Command line argument parsing using clap
    let matches = Command::new("Rust BAM Processor")
        .version("0.1")
        .author("Your Name <your.email@example.com>")
        .about("Processes BAM files")
        .arg(
            Arg::new("vcf-file")
                .long("vcf-file")
                .value_parser(clap::value_parser!(String))
                .help("VCF file to process")
                .required(true),
        )
        .arg(
            Arg::new("bam-files")
                .long("bam-files")
                .value_parser(clap::value_parser!(String))
                .help("Comma-separated list of BAM files to process")
                .required(true),
        )
        .get_matches();

    // Extract arguments
    let vcf_file = matches.get_one::<String>("vcf-file").unwrap();
    let bam_files = matches.get_one::<String>("bam-files").unwrap();

    // Process VCF file
    if let Err(e) = read_vcf_file(vcf_file) {
        error!("Error reading VCF file: {}", e);
    }

    // Process BAM files
    let bam_file_paths: Vec<&str> = bam_files.split(',').collect();
    
    println!("Processing {} BAM files...", bam_file_paths.len());
    
    bam_file_paths.par_iter().for_each(|&bam_file| {
        if let Err(e) = read_bam_file(bam_file) {
            error!("Error reading BAM file {}: {}", bam_file, e);
        }
    });
    
    println!("Processing completed successfully");
    Ok(())
}