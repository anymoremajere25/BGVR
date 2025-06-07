// src/bin/tss_calculator.rs
use clap::Parser;
use serde_json::json;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

#[derive(Parser)]
#[command(name = "tss_calculator")]
#[command(about = "Calculate TSS enrichment scores for ATAC-Seq data")]
pub struct Args {
    #[arg(long, help = "Input BAM file")]
    bam: String,
    
    #[arg(long, help = "TSS sites BED file")]
    tss_bed: String,
    
    #[arg(long, default_value = "2000", help = "Window size around TSS")]
    window_size: i32,
    
    #[arg(long, help = "Output JSON file")]
    output_json: String,
    
    #[arg(long, help = "Output profile BED file")]
    output_profile: String,
    
    #[arg(long, help = "Output plot PNG file")]
    output_plot: String,
    
    #[arg(long, help = "Sample name")]
    sample_name: String,
    
    #[arg(long, default_value = "2", help = "Number of threads")]
    threads: usize,
    
    #[arg(long, help = "Log file path")]
    log_message(&mut log_writer, "Starting TSS enrichment calculation")?;
    log_message(&mut log_writer, &format!("BAM file: {}", args.bam))?;
    log_message(&mut log_writer, &format!("TSS BED file: {}", args.tss_bed))?;
    log_message(&mut log_writer, &format!("Window size: {}", args.window_size))?;
    log_message(&mut log_writer, &format!("Sample: {}", args.sample_name))?;
    
    // Step 1: Load TSS sites from BED file
    log_message(&mut log_writer, "Loading TSS sites from BED file...")?;
    let tss_sites = load_tss_sites(&args.tss_bed, &mut log_writer)?;
    log_message(&mut log_writer, &format!("Loaded {} TSS sites", tss_sites.len()))?;
    
    // Step 2: Calculate coverage around TSS sites
    log_message(&mut log_writer, "Calculating coverage around TSS sites...")?;
    let coverage_data = calculate_tss_coverage(&args.bam, &tss_sites, args.window_size, &mut log_writer)?;
    
    // Step 3: Compute enrichment metrics
    log_message(&mut log_writer, "Computing enrichment metrics...")?;
    let enrichment_result = compute_enrichment_metrics(coverage_data, args.window_size, tss_sites.len())?;
    
    log_message(&mut log_writer, &format!("TSS enrichment score: {:.2}", enrichment_result.enrichment_score))?;
    log_message(&mut log_writer, &format!("Signal-to-noise ratio: {:.2}", enrichment_result.signal_to_noise))?;
    log_message(&mut log_writer, &format!("Quality grade: {}", enrichment_result.get_quality_grade()))?;
    
    // Step 4: Write outputs
    log_message(&mut log_writer, "Writing output files...")?;
    
    // Write JSON results
    let json_output = enrichment_result.to_json(&args.sample_name);
    std::fs::write(&args.output_json, serde_json::to_string_pretty(&json_output)?)?;
    log_message(&mut log_writer, &format!("JSON results written to: {}", args.output_json))?;
    
    // Write coverage profile BED
    write_coverage_profile(&enrichment_result, &args.output_profile, &mut log_writer)?;
    
    // Generate plot
    generate_tss_plot(&enrichment_result, &args.output_plot, &args.sample_name, &mut log_writer)?;
    
    log_message(&mut log_writer, "TSS enrichment calculation completed successfully")?;
    
    Ok(())
}

fn log_message(writer: &mut Option<BufWriter<File>>, message: &str) -> std::io::Result<()> {
    let timestamp = chrono::Utc::now().format("%Y-%m-%d %H:%M:%S");
    let log_line = format!("[{}] {}", timestamp, message);
    
    println!("{}", log_line);
    
    if let Some(w) = writer {
        writeln!(w, "{}", log_line)?;
        w.flush()?;
    }
    
    Ok(())
}

fn load_tss_sites(bed_file: &str, log_writer: &mut Option<BufWriter<File>>) -> Result<Vec<TSSsite>, Box<dyn std::error::Error>> {
    let file = File::open(bed_file)?;
    let reader = BufReader::new(file);
    let mut tss_sites = Vec::new();
    
    for (line_num, line) in reader.lines().enumerate() {
        let line = line?;
        
        // Skip header lines and comments
        if line.starts_with('#') || line.starts_with("track") || line.trim().is_empty() {
            continue;
        }
        
        let fields: Vec<&str> = line.split('\t').collect();
        
        if fields.len() < 3 {
            log_message(log_writer, &format!("Warning: Invalid BED line {}: {}", line_num + 1, line))?;
            continue;
        }
        
        let chromosome = fields[0].to_string();
        let start: i32 = fields[1].parse().map_err(|_| {
            format!("Invalid start coordinate at line {}: {}", line_num + 1, fields[1])
        })?;
        let end: i32 = fields[2].parse().map_err(|_| {
            format!("Invalid end coordinate at line {}: {}", line_num + 1, fields[2])
        })?;
        
        let name = if fields.len() > 3 {
            fields[3].to_string()
        } else {
            format!("TSS_{}", line_num + 1)
        };
        
        let score = if fields.len() > 4 {
            fields[4].parse().unwrap_or(0.0)
        } else {
            0.0
        };
        
        let strand = if fields.len() > 5 {
            fields[5].chars().next().unwrap_or('+')
        } else {
            '+'
        };
        
        tss_sites.push(TSSsite {
            chromosome,
            start,
            end,
            name,
            score,
            strand,
        });
    }
    
    Ok(tss_sites)
}

fn calculate_tss_coverage(bam_file: &str, tss_sites: &[TSSsite], window_size: i32, log_writer: &mut Option<BufWriter<File>>) -> Result<Vec<f64>, Box<dyn std::error::Error>> {
    let half_window = window_size / 2;
    let mut coverage_sum = vec![0.0; window_size as usize];
    let mut processed_sites = 0;
    
    for (i, tss) in tss_sites.iter().enumerate() {
        if i % 1000 == 0 {
            log_message(log_writer, &format!("Processing TSS site {} of {}", i + 1, tss_sites.len()))?;
        }
        
        // Calculate TSS center
        let tss_center = (tss.start + tss.end) / 2;
        let region_start = tss_center - half_window;
        let region_end = tss_center + half_window;
        
        // Get coverage for this region using samtools
        let coverage = get_region_coverage(bam_file, &tss.chromosome, region_start, region_end)?;
        
        if coverage.len() == window_size as usize {
            for (j, &cov) in coverage.iter().enumerate() {
                coverage_sum[j] += cov;
            }
            processed_sites += 1;
        }
    }
    
    // Average the coverage across all sites
    if processed_sites > 0 {
        for cov in &mut coverage_sum {
            *cov /= processed_sites as f64;
        }
    }
    
    log_message(log_writer, &format!("Successfully processed {} TSS sites", processed_sites))?;
    
    Ok(coverage_sum)
}

fn get_region_coverage(bam_file: &str, chromosome: &str, start: i32, end: i32) -> Result<Vec<f64>, Box<dyn std::error::Error>> {
    use std::process::{Command, Stdio};
    
    let region = format!("{}:{}-{}", chromosome, start, end);
    
    let mut cmd = Command::new("samtools");
    cmd.args(&[
        "depth",
        "-r", &region,
        bam_file,
    ]);
    
    cmd.stdout(Stdio::piped());
    let output = cmd.output()?;
    
    if !output.status.success() {
        return Err(format!("samtools depth failed for region {}", region).into());
    }
    
    let depth_output = String::from_utf8(output.stdout)?;
    let mut coverage = vec![0.0; (end - start) as usize];
    
    for line in depth_output.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 3 {
            if let (Ok(pos), Ok(depth)) = (fields[1].parse::<i32>(), fields[2].parse::<f64>()) {
                let index = (pos - start) as usize;
                if index < coverage.len() {
                    coverage[index] = depth;
                }
            }
        }
    }
    
    Ok(coverage)
}

fn compute_enrichment_metrics(coverage_profile: Vec<f64>, window_size: i32, total_sites: usize) -> Result<TSSEnrichmentResult, Box<dyn std::error::Error>> {
    let half_window = window_size / 2;
    let center_region_size = 100; // ±50bp around TSS
    let center_start = (half_window - center_region_size / 2) as usize;
    let center_end = (half_window + center_region_size / 2) as usize;
    
    // Calculate max signal in center region
    let max_signal = coverage_profile[center_start..center_end]
        .iter()
        .fold(0.0, |max, &val| max.max(val));
    
    // Calculate background signal (flanking regions)
    let flank_size = 500; // 500bp flanks
    let left_flank_end = flank_size;
    let right_flank_start = coverage_profile.len() - flank_size;
    
    let left_background: f64 = coverage_profile[0..left_flank_end].iter().sum::<f64>() / flank_size as f64;
    let right_background: f64 = coverage_profile[right_flank_start..].iter().sum::<f64>() / flank_size as f64;
    let background_signal = (left_background + right_background) / 2.0;
    
    // Calculate enrichment score
    let enrichment_score = if background_signal > 0.0 {
        max_signal / background_signal
    } else {
        max_signal
    };
    
    // Calculate signal-to-noise ratio
    let center_mean: f64 = coverage_profile[center_start..center_end].iter().sum::<f64>() 
        / (center_end - center_start) as f64;
    let signal_to_noise = if background_signal > 0.0 {
        center_mean / background_signal
    } else {
        center_mean
    };
    
    Ok(TSSEnrichmentResult {
        enrichment_score,
        signal_to_noise,
        max_signal,
        background_signal,
        total_sites,
        coverage_profile,
        window_size,
    })
}

fn write_coverage_profile(result: &TSSEnrichmentResult, output_file: &str, log_writer: &mut Option<BufWriter<File>>) -> Result<(), Box<dyn std::error::Error>> {
    let mut writer = BufWriter::new(File::create(output_file)?);
    let half_window = result.window_size / 2;
    
    writeln!(writer, "# TSS Coverage Profile")?;
    writeln!(writer, "# Position_relative_to_TSS\tCoverage")?;
    
    for (i, &coverage) in result.coverage_profile.iter().enumerate() {
        let position = i as i32 - half_window;
        writeln!(writer, "{}\t{:.6}", position, coverage)?;
    }
    
    writer.flush()?;
    log_message(log_writer, &format!("Coverage profile written to: {}", output_file))?;
    
    Ok(())
}

fn generate_tss_plot(result: &TSSEnrichmentResult, output_file: &str, sample_name: &str, log_writer: &mut Option<BufWriter<File>>) -> Result<(), Box<dyn std::error::Error>> {
    // Create a simple Python script to generate the plot
    let python_script = format!(r#"
import matplotlib.pyplot as plt
import numpy as np

# Data
positions = np.arange({}, {})
coverage = {}

# Create plot
plt.figure(figsize=(10, 6))
plt.plot(positions, coverage, linewidth=2, color='blue')
plt.axvline(x=0, color='red', linestyle='--', alpha=0.7, label='TSS')
plt.xlabel('Distance from TSS (bp)')
plt.ylabel('Average Coverage')
plt.title('TSS Enrichment Profile - {}\nEnrichment Score: {:.2f}')
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig('{}', dpi=300, bbox_inches='tight')
plt.close()
print("Plot saved successfully")
"#, 
        -result.window_size / 2,
        result.window_size / 2,
        format!("{:?}", result.coverage_profile),
        sample_name,
        result.enrichment_score,
        output_file
    );
    
    // Write and execute Python script
    let script_file = format!("{}.py", output_file.trim_end_matches(".png"));
    std::fs::write(&script_file, python_script)?;
    
    let output = std::process::Command::new("python3")
        .arg(&script_file)
        .output()?;
    
    if output.status.success() {
        std::fs::remove_file(&script_file)?; // Clean up script
        log_message(log_writer, &format!("TSS plot generated: {}", output_file))?;
    } else {
        let error_msg = String::from_utf8_lossy(&output.stderr);
        log_message(log_writer, &format!("Plot generation warning: {}", error_msg))?;
        
        // Create a simple text-based plot as fallback
        create_text_plot(result, output_file, sample_name)?;
    }
    
    Ok(())
}

fn create_text_plot(result: &TSSEnrichmentResult, output_file: &str, sample_name: &str) -> Result<(), Box<dyn std::error::Error>> {
    let mut writer = BufWriter::new(File::create(output_file.replace(".png", "_plot.txt"))?);
    
    writeln!(writer, "TSS Enrichment Profile - {}", sample_name)?;
    writeln!(writer, "Enrichment Score: {:.2}", result.enrichment_score)?;
    writeln!(writer, "Signal-to-Noise: {:.2}", result.signal_to_noise)?;
    writeln!(writer, "=")?;
    
    let max_coverage = result.coverage_profile.iter().fold(0.0, |max, &val| max.max(val));
    let scale_factor = 50.0 / max_coverage; // Scale to 50 characters width
    
    for (i, &coverage) in result.coverage_profile.iter().enumerate() {
        let position = i as i32 - result.window_size / 2;
        let bar_length = (coverage * scale_factor) as usize;
        let bar = "#".repeat(bar_length);
        
        if position % 200 == 0 { // Show every 200bp
            writeln!(writer, "{:6}: {} ({:.2})", position, bar, coverage)?;
        }
    }
    
    writer.flush()?;
    Ok(())
}file: Option<String>,
}

#[derive(Debug, Clone)]
pub struct TSSsite {
    chromosome: String,
    start: i32,
    end: i32,
    name: String,
    score: f64,
    strand: char,
}

#[derive(Debug)]
pub struct TSSEnrichmentResult {
    enrichment_score: f64,
    signal_to_noise: f64,
    max_signal: f64,
    background_signal: f64,
    total_sites: usize,
    coverage_profile: Vec<f64>,
    window_size: i32,
}

impl TSSEnrichmentResult {
    fn to_json(&self, sample_name: &str) -> serde_json::Value {
        json!({
            "sample": sample_name,
            "tss_enrichment": self.enrichment_score,
            "signal_to_noise": self.signal_to_noise,
            "max_signal": self.max_signal,
            "background_signal": self.background_signal,
            "total_tss_sites": self.total_sites,
            "window_size": self.window_size,
            "quality_grade": self.get_quality_grade(),
            "analysis_timestamp": chrono::Utc::now().to_rfc3339()
        })
    }
    
    fn get_quality_grade(&self) -> &str {
        match self.enrichment_score {
            x if x >= 10.0 => "Excellent",
            x if x >= 7.0 => "Good",
            x if x >= 5.0 => "Acceptable",
            x if x >= 3.0 => "Poor",
            _ => "Failed"
        }
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    
    // Initialize logging
    let mut log_writer = if let Some(log_file) = &args.log_file {
        Some(BufWriter::new(File::create(log_file)?))
    } else {
        None
    };
    
    log_