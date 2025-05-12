use statrs::distribution::{ChiSquared, ContinuousCDF};
use polars::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
// 'rust-htslib' for VCF/BCF I/O
// 'statrs' for statistical distributions
// 'polars' for data frame operations

fn chi_square_hw(aa: f64, ab: f64, bb: f64, p: f64) -> f64 {
    let total = aa + ab + bb;
    if total == 0.0 {
        return 1.0;
    }
    let q = 1.0 - p;
    let expected_aa = p * p * total;
    let expected_ab = 2.0 * p * q * total;
    let expected_bb = q * q * total;
    let chi_sq = (aa - expected_aa).powi(2) / expected_aa
               + (ab - expected_ab).powi(2) / expected_ab
               + (bb - expected_bb).powi(2) / expected_bb;
    let dist = ChiSquared::new(1.0).unwrap();
    1.0 - dist.cdf(chi_sq)
}

// Custom parser for VCF files with potential formatting issues
fn process_vcf_file(vcf_path: &Path, start_pos: u64, end_pos: u64) -> Result<DataFrame, Box<dyn std::error::Error + Send + Sync>> {
    // Open the file
    let file = File::open(vcf_path)?;
    let reader = BufReader::new(file);
    
    // Parse the file line by line
    let mut records_data: Vec<(String, u64, String, String, f64)> = Vec::new();
    let mut sample_indices = Vec::new();
    
    for line in reader.lines() {
        let line = line?;
        
        // Skip comment lines
        if line.starts_with("##") {
            continue;
        }
        
        // Process header line
        if line.starts_with("#CHROM") {
            let columns: Vec<&str> = line.split_whitespace().collect();
            
            // Find sample columns (FORMAT column is at index 8, samples start at 9)
            if columns.len() > 9 {
                sample_indices = (9..columns.len()).collect();
            }
            continue;
        }
        
        // Process data lines
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 10 {
            // Skip malformed lines
            continue;
        }
        
        // Extract basic variant information
        let chrom = fields[0].to_string();
        let pos: u64 = fields[1].parse()?;
        
        // Apply position filter
        if pos < start_pos || pos > end_pos {
            continue;
        }
        
        let ref_allele = fields[3].to_string();
        let alt_allele = fields[4].to_string();
        
        // Find FORMAT field index (typically 8)
        let format_field = fields[8];
        let format_columns: Vec<&str> = format_field.split(':').collect();
        let gt_index = format_columns.iter().position(|&x| x == "GT");
        
        if gt_index.is_none() {
            // Skip if no GT field
            continue;
        }
        
        // Count genotypes
        let mut count_aa = 0.0;
        let mut count_ab = 0.0;
        let mut count_bb = 0.0;
        
        for &sample_idx in &sample_indices {
            if sample_idx >= fields.len() {
                continue;
            }
            
            let sample_data = fields[sample_idx];
            let sample_fields: Vec<&str> = sample_data.split(':').collect();
            
            let gt_idx = gt_index.unwrap();
            if gt_idx >= sample_fields.len() {
                continue;
            }
            
            let genotype = sample_fields[gt_idx];
            match genotype {
                "0/0" | "0|0" => count_aa += 1.0,
                "0/1" | "1/0" | "0|1" | "1|0" => count_ab += 1.0,
                "1/1" | "1|1" => count_bb += 1.0,
                _ => {} // Skip other genotypes like ./. or multi-allelic
            }
        }
        
        let total = count_aa + count_ab + count_bb;
        if total == 0.0 {
            continue;
        }
        
        // Calculate allele frequency and HW equilibrium
        let p = ((count_aa * 2.0) + count_ab) / (2.0 * total);
        let hw_p = chi_square_hw(count_aa, count_ab, count_bb, p);
        
        records_data.push((
            chrom,
            pos,
            ref_allele,
            alt_allele,
            hw_p,
        ));
    }
    
    // Create DataFrame from the parsed data
    let chrom_series = Series::new("chrom".into(), records_data.iter().map(|x| x.0.clone()).collect::<Vec<String>>());
    let pos_series = Series::new("pos".into(), records_data.iter().map(|x| x.1).collect::<Vec<u64>>());
    let ref_series = Series::new("ref_allele".into(), records_data.iter().map(|x| x.2.clone()).collect::<Vec<String>>());
    let alt_series = Series::new("alt_allele".into(), records_data.iter().map(|x| x.3.clone()).collect::<Vec<String>>());
    let hw_series = Series::new("hw_pvalue".into(), records_data.iter().map(|x| x.4).collect::<Vec<f64>>());
    
    let df = DataFrame::new(vec![
        chrom_series.into(),
        pos_series.into(),
        ref_series.into(),
        alt_series.into(),
        hw_series.into(),
    ])?;
    
    Ok(df)
}

fn main() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} <vcf_file> [start_pos] [end_pos]", args[0]);
        std::process::exit(1);
    }
    
    let vcf_path = Path::new(&args[1]);
    let start_pos: u64 = if args.len() > 2 { args[2].parse()? } else { 0 };
    let end_pos: u64 = if args.len() > 3 { args[3].parse()? } else { u64::MAX };
    
    println!("Processing VCF file: {}", vcf_path.display());
    println!("Position range: {} - {}", start_pos, end_pos);
    
    // Process the VCF file using our custom parser
    match process_vcf_file(vcf_path, start_pos, end_pos) {
        Ok(mut df) => {  // Make df mutable
            println!("Analysis complete. Results:");
            println!("{}", df);
            
            // Save results to CSV using the proper method with mutable reference
            let output_path = format!("{}.hw_results.csv", vcf_path.display());
            match CsvWriter::new(File::create(&output_path)?)
                .finish(&mut df) {  // Pass a mutable reference
                Ok(_) => println!("Results saved to {}", output_path),
                Err(e) => eprintln!("Failed to save results: {}", e),
            }
        },
        Err(e) => {
            eprintln!("Error processing VCF file: {}", e);
            return Err(e);
        }
    }
    
    Ok(())
}