// src/bin/atac_aligner.rs
use clap::Parser;
use serde_json::json;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::process::{Command, Stdio};
use std::time::Instant;

#[derive(Parser)]
#[command(name = "atac_aligner")]
#[command(about = "ATAC-Seq read aligner with quality filtering")]
pub struct Args {
    #[arg(long, help = "Reference genome FASTA file")]
    reference: String,
    
    #[arg(long, help = "Read 1 FASTQ file")]
    read1: String,
    
    #[arg(long, help = "Read 2 FASTQ file")]
    read2: String,
    
    #[arg(long, help = "Output BAM file")]
    output: String,
    
    #[arg(long, default_value = "4", help = "Number of threads")]
    threads: usize,
    
    #[arg(long, default_value = "30", help = "Minimum mapping quality")]
    quality_threshold: u8,
    
    #[arg(long, help = "Remove duplicate reads")]
    remove_duplicates: bool,
    
    #[arg(long, help = "Remove mitochondrial reads")]
    remove_mitochondrial: bool,
    
    #[arg(long, help = "Log file path")]
    log_file: Option<String>,
    
    #[arg(long, help = "Statistics JSON output")]
    stats_json: Option<String>,
}

#[derive(Debug)]
pub struct AlignmentStats {
    total_reads: u64,
    mapped_reads: u64,
    high_quality_reads: u64,
    duplicate_reads: u64,
    mitochondrial_reads: u64,
    final_reads: u64,
    alignment_time: f64,
}

impl AlignmentStats {
    fn new() -> Self {
        Self {
            total_reads: 0,
            mapped_reads: 0,
            high_quality_reads: 0,
            duplicate_reads: 0,
            mitochondrial_reads: 0,
            final_reads: 0,
            alignment_time: 0.0,
        }
    }
    
    fn to_json(&self, sample_name: &str) -> serde_json::Value {
        json!({
            "sample": sample_name,
            "total_reads": self.total_reads,
            "mapped_reads": self.mapped_reads,
            "high_quality_reads": self.high_quality_reads,
            "duplicate_reads": self.duplicate_reads,
            "mitochondrial_reads": self.mitochondrial_reads,
            "final_reads": self.final_reads,
            "mapping_rate": if self.total_reads > 0 { 
                (self.mapped_reads as f64 / self.total_reads as f64) * 100.0 
            } else { 0.0 },
            "quality_rate": if self.mapped_reads > 0 { 
                (self.high_quality_reads as f64 / self.mapped_reads as f64) * 100.0 
            } else { 0.0 },
            "duplicate_rate": if self.mapped_reads > 0 { 
                (self.duplicate_reads as f64 / self.mapped_reads as f64) * 100.0 
            } else { 0.0 },
            "alignment_time_seconds": self.alignment_time
        })
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let start_time = Instant::now();
    
    // Initialize logging
    let mut log_writer = if let Some(log_file) = &args.log_file {
        Some(BufWriter::new(File::create(log_file)?))
    } else {
        None
    };
    
    log_message(&mut log_writer, "Starting ATAC-Seq alignment pipeline")?;
    log_message(&mut log_writer, &format!("Reference: {}", args.reference))?;
    log_message(&mut log_writer, &format!("Read1: {}", args.read1))?;
    log_message(&mut log_writer, &format!("Read2: {}", args.read2))?;
    log_message(&mut log_writer, &format!("Output: {}", args.output))?;
    log_message(&mut log_writer, &format!("Threads: {}", args.threads))?;
    
    let mut stats = AlignmentStats::new();
    
    // Step 1: Count input reads
    log_message(&mut log_writer, "Counting input reads...")?;
    stats.total_reads = count_fastq_reads(&args.read1)?;
    log_message(&mut log_writer, &format!("Total read pairs: {}", stats.total_reads))?;
    
    // Step 2: Align reads using BWA-MEM2 (or similar aligner)
    log_message(&mut log_writer, "Aligning reads with BWA-MEM2...")?;
    let sam_file = format!("{}.sam", args.output.trim_end_matches(".bam"));
    align_reads(&args, &sam_file, &mut log_writer)?;
    
    // Step 3: Convert SAM to BAM and filter
    log_message(&mut log_writer, "Converting to BAM and filtering...")?;
    let filtered_bam = format!("{}.filtered.bam", args.output.trim_end_matches(".bam"));
    convert_and_filter(&sam_file, &filtered_bam, &args, &mut stats, &mut log_writer)?;
    
    // Step 4: Remove duplicates if requested
    let dedup_bam = if args.remove_duplicates {
        log_message(&mut log_writer, "Removing duplicates...")?;
        let dedup_file = format!("{}.dedup.bam", args.output.trim_end_matches(".bam"));
        remove_duplicates(&filtered_bam, &dedup_file, &mut stats, &mut log_writer)?;
        dedup_file
    } else {
        filtered_bam
    };
    
    // Step 5: Remove mitochondrial reads if requested
    let final_bam = if args.remove_mitochondrial {
        log_message(&mut log_writer, "Removing mitochondrial reads...")?;
        remove_mitochondrial_reads(&dedup_bam, &args.output, &mut stats, &mut log_writer)?;
        args.output.clone()
    } else {
        // Just rename the file to final output
        std::fs::rename(&dedup_bam, &args.output)?;
        args.output.clone()
    };
    
    // Step 6: Sort and index final BAM
    log_message(&mut log_writer, "Sorting final BAM file...")?;
    sort_bam(&final_bam, &mut log_writer)?;
    
    // Calculate final statistics
    stats.alignment_time = start_time.elapsed().as_secs_f64();
    stats.final_reads = count_bam_reads(&args.output)?;
    
    log_message(&mut log_writer, &format!("Final read count: {}", stats.final_reads))?;
    log_message(&mut log_writer, &format!("Alignment completed in {:.2} seconds", stats.alignment_time))?;
    
    // Write statistics JSON
    if let Some(stats_file) = &args.stats_json {
        let sample_name = Path::new(&args.output)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown");
        let stats_json = stats.to_json(sample_name);
        std::fs::write(stats_file, serde_json::to_string_pretty(&stats_json)?)?;
        log_message(&mut log_writer, &format!("Statistics written to: {}", stats_file))?;
    }
    
    // Cleanup temporary files
    cleanup_temp_files(&sam_file, &mut log_writer)?;
    
    log_message(&mut log_writer, "ATAC-Seq alignment pipeline completed successfully")?;
    
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

fn count_fastq_reads(fastq_file: &str) -> Result<u64, Box<dyn std::error::Error>> {
    let output = if fastq_file.ends_with(".gz") {
        Command::new("zcat")
            .arg(fastq_file)
            .stdout(Stdio::piped())
            .spawn()?
            .wait_with_output()?
    } else {
        Command::new("cat")
            .arg(fastq_file)
            .stdout(Stdio::piped())
            .spawn()?
            .wait_with_output()?
    };
    
    let line_count = String::from_utf8(output.stdout)?
        .lines()
        .count() as u64;
    
    Ok(line_count / 4) // FASTQ has 4 lines per read
}

fn align_reads(args: &Args, output_sam: &str, log_writer: &mut Option<BufWriter<File>>) -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::new("bwa-mem2");
    cmd.args(&[
        "mem",
        "-t", &args.threads.to_string(),
        "-M", // Mark shorter split hits as secondary
        &args.reference,
        &args.read1,
        &args.read2,
    ]);
    
    cmd.stdout(File::create(output_sam)?);
    
    let output = cmd.output()?;
    
    if !output.status.success() {
        let error_msg = String::from_utf8_lossy(&output.stderr);
        log_message(log_writer, &format!("BWA-MEM2 error: {}", error_msg))?;
        return Err(format!("BWA-MEM2 failed with exit code: {}", output.status).into());
    }
    
    Ok(())
}

fn convert_and_filter(sam_file: &str, output_bam: &str, args: &Args, stats: &mut AlignmentStats, log_writer: &mut Option<BufWriter<File>>) -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::new("samtools");
    cmd.args(&[
        "view",
        "-b", // Output BAM
        "-h", // Include header
        "-F", "4", // Exclude unmapped reads
        "-q", &args.quality_threshold.to_string(), // Minimum mapping quality
        "-@", &args.threads.to_string(),
        sam_file,
        "-o", output_bam,
    ]);
    
    let output = cmd.output()?;
    
    if !output.status.success() {
        let error_msg = String::from_utf8_lossy(&output.stderr);
        log_message(log_writer, &format!("Samtools filter error: {}", error_msg))?;
        return Err("Samtools filtering failed".into());
    }
    
    // Update statistics
    stats.mapped_reads = count_bam_reads(output_bam)?;
    stats.high_quality_reads = stats.mapped_reads;
    
    log_message(log_writer, &format!("Mapped reads after filtering: {}", stats.mapped_reads))?;
    
    Ok(())
}

fn remove_duplicates(input_bam: &str, output_bam: &str, stats: &mut AlignmentStats, log_writer: &mut Option<BufWriter<File>>) -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::new("samtools");
    cmd.args(&[
        "rmdup",
        input_bam,
        output_bam,
    ]);
    
    let output = cmd.output()?;
    
    if !output.status.success() {
        let error_msg = String::from_utf8_lossy(&output.stderr);
        log_message(log_writer, &format!("Duplicate removal error: {}", error_msg))?;
        return Err("Duplicate removal failed".into());
    }
    
    let reads_after_dedup = count_bam_reads(output_bam)?;
    stats.duplicate_reads = stats.high_quality_reads - reads_after_dedup;
    
    log_message(log_writer, &format!("Duplicates removed: {}", stats.duplicate_reads))?;
    
    Ok(())
}

fn remove_mitochondrial_reads(input_bam: &str, output_bam: &str, stats: &mut AlignmentStats, log_writer: &mut Option<BufWriter<File>>) -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::new("samtools");
    cmd.args(&[
        "view",
        "-b",
        "-h",
        input_bam,
        // Exclude common mitochondrial chromosome names
        "-L", "/dev/stdin",
    ]);
    
    // Create a temporary BED file excluding mitochondrial regions
    let exclude_regions = "chrM\t0\t20000\nchrMT\t0\t20000\nMT\t0\t20000\nM\t0\t20000\n";
    cmd.stdin(Stdio::piped());
    
    let mut child = cmd.spawn()?;
    
    if let Some(stdin) = child.stdin.as_mut() {
        stdin.write_all(exclude_regions.as_bytes())?;
    }
    
    let output = child.wait_with_output()?;
    
    if output.status.success() {
        std::fs::write(output_bam, &output.stdout)?;
        let reads_after_mito_filter = count_bam_reads(output_bam)?;
        stats.mitochondrial_reads = count_bam_reads(input_bam)? - reads_after_mito_filter;
        
        log_message(log_writer, &format!("Mitochondrial reads removed: {}", stats.mitochondrial_reads))?;
    } else {
        return Err("Mitochondrial read removal failed".into());
    }
    
    Ok(())
}

fn sort_bam(bam_file: &str, log_writer: &mut Option<BufWriter<File>>) -> Result<(), Box<dyn std::error::Error>> {
    let temp_sorted = format!("{}.temp_sorted", bam_file);
    
    let mut cmd = Command::new("samtools");
    cmd.args(&[
        "sort",
        "-o", &temp_sorted,
        bam_file,
    ]);
    
    let output = cmd.output()?;
    
    if output.status.success() {
        std::fs::rename(&temp_sorted, bam_file)?;
        log_message(log_writer, "BAM file sorted successfully")?;
    } else {
        let error_msg = String::from_utf8_lossy(&output.stderr);
        log_message(log_writer, &format!("BAM sorting error: {}", error_msg))?;
        return Err("BAM sorting failed".into());
    }
    
    Ok(())
}

fn count_bam_reads(bam_file: &str) -> Result<u64, Box<dyn std::error::Error>> {
    let output = Command::new("samtools")
        .args(&["view", "-c", bam_file])
        .output()?;
    
    if output.status.success() {
        let count_str = String::from_utf8(output.stdout)?;
        Ok(count_str.trim().parse()?)
    } else {
        Err("Failed to count BAM reads".into())
    }
}

fn cleanup_temp_files(sam_file: &str, log_writer: &mut Option<BufWriter<File>>) -> Result<(), Box<dyn std::error::Error>> {
    if Path::new(sam_file).exists() {
        std::fs::remove_file(sam_file)?;
        log_message(log_writer, "Temporary SAM file removed")?;
    }
    Ok(())
}