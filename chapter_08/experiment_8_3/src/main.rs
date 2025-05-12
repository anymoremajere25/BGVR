use anyhow::{anyhow, Context, Result};
use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use std::{
    collections::HashMap,
    fs::File,
    io::{BufReader, BufWriter},
    path::{Path, PathBuf},
    sync::{Arc, Mutex},
    time::Instant,
};
use thiserror::Error;
use tracing::{debug, info, warn, Level};
use polars::prelude::*;
use serde::{Deserialize, Serialize};
use std::str::FromStr;

/// Custom error type
#[derive(Error, Debug)]
pub enum CallerError {
    #[error("No variants detected")]
    NoVariants,
    #[error("I/O error: {0}")]
    IoError(#[from] std::io::Error),
    #[error("{0}")]
    Other(String),
}

/// Command-line arguments
#[derive(Parser, Debug)]
#[clap(author, version, about = "Fast variant caller implemented in Rust")]
struct Cli {
    /// BAM file path
    #[arg(short, long)]
    bam: PathBuf,

    /// Reference FASTA (index .fai required)
    #[arg(short, long)]
    fasta: Option<PathBuf>,

    /// Region filter (chr:start-end)
    #[arg(short, long)]
    region: Option<String>,

    /// Output Parquet
    #[arg(short, long)]
    out: PathBuf,

    /// Minimum read depth
    #[arg(long, default_value_t = 8)]
    min_depth: usize,

    /// Minimum genotype quality
    #[arg(long, default_value_t = 20.0)]
    min_gq: f32,

    /// Minimum mapping quality
    #[arg(long, default_value_t = 20)]
    min_mapq: u8,

    /// Minimum base quality
    #[arg(long, default_value_t = 20)]
    min_baseq: u8,

    /// Threads (0=auto)
    #[arg(short, long, default_value_t = 0)]
    threads: usize,

    /// Verbose logging
    #[arg(short, long)]
    verbose: bool,

    /// Export stats JSON
    #[arg(long)]
    stats: Option<PathBuf>,
}

/// Pileup entry for a single position
#[derive(Debug, Default)]
struct PileupEntry {
    depth: u32,
    base_counts: HashMap<char, u32>,
    forward_strands: HashMap<char, u32>,
    reverse_strands: HashMap<char, u32>,
    total_mapq: u32,
    total_baseq: u32,
}

impl PileupEntry {
    fn add_base(&mut self, base: char, is_forward: bool, mapq: u8, baseq: u8) {
        self.depth += 1;
        *self.base_counts.entry(base).or_insert(0) += 1;
        if is_forward { *self.forward_strands.entry(base).or_insert(0) += 1; }
        else { *self.reverse_strands.entry(base).or_insert(0) += 1; }
        self.total_mapq += mapq as u32;
        self.total_baseq += baseq as u32;
    }
    fn get_calls(&self, ref_base: char, min_depth: usize, min_gq: f32, min_mapq: u8, min_baseq: u8) -> Vec<Call> {
        let mut calls = Vec::new();
        if self.depth < min_depth as u32 { return calls; }
        let mapq_avg = self.total_mapq as f32 / self.depth as f32;
        let baseq_avg = self.total_baseq as f32 / self.depth as f32;
        if mapq_avg < min_mapq as f32 || baseq_avg < min_baseq as f32 { return calls; }
        let ref_count = *self.base_counts.get(&ref_base).unwrap_or(&0);
        for (&alt, &count) in &self.base_counts {
            if alt == ref_base || count == 0 { continue; }
            let vaf = count as f32 / self.depth as f32;
            let gq = -10.0 * (0.5 - (vaf - 0.5).abs()).log10();
            if gq < min_gq { continue; }
            let fwd_alt = *self.forward_strands.get(&alt).unwrap_or(&0) as f32;
            let rev_alt = *self.reverse_strands.get(&alt).unwrap_or(&0) as f32;
            let fwd_ref = *self.forward_strands.get(&ref_base).unwrap_or(&0) as f32;
            let rev_ref = *self.reverse_strands.get(&ref_base).unwrap_or(&0) as f32;
            let strand_bias = if fwd_alt + rev_alt > 0.0 && fwd_ref + rev_ref > 0.0 {
                let diff = (fwd_alt/(fwd_alt+rev_alt) - fwd_ref/(fwd_ref+rev_ref)).abs();
                1.0 - diff
            } else { 0.0 };
            calls.push(Call {
                chrom: String::new(), pos: 0, ref_base, alt_base: alt,
                depth: self.depth, ref_count, alt_count: count,
                gq, mapq_avg, baseq_avg, vaf, strand_bias,
            });
        }
        calls
    }
}

/// Simplified BAM reader stub
struct SimpleBamReader { _hdr: SimpleHeader }
impl SimpleBamReader {
    fn new(_p: &Path) -> Result<Self> { Ok(SimpleBamReader{ _hdr: SimpleHeader::default() }) }
    fn read_header(&mut self) -> Result<SimpleHeader> { Ok(self._hdr.clone()) }
}
#[derive(Clone, Default)] struct SimpleHeader;
impl SimpleHeader { fn reference_sequences(&self) -> HashMap<String,usize> {
    let mut m = HashMap::new(); m.insert("chr1".into(),248_956_422); m
}}

/// Region struct
#[derive(Debug, Clone)] struct Region { name:String, start:Option<usize>, end:Option<usize> }
impl FromStr for Region {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let parts = s.split(':').collect::<Vec<_>>();
        let name = parts[0].to_string();

        if parts.len() == 1 {
            return Ok(Region { name, start: None, end: None });
        }

        let range = parts[1];
        let bounds = range.split('-').collect::<Vec<_>>();
        let start = bounds.get(0).and_then(|b| b.parse::<usize>().ok());
        let end = bounds.get(1).and_then(|b| b.parse::<usize>().ok());

        Ok(Region { name, start, end })
    }
}

/// Variant call record
#[derive(Serialize, Deserialize, Debug, Clone)] struct Call {
    chrom:String,pos:i64,ref_base:char,alt_base:char,
    depth:u32,ref_count:u32,alt_count:u32,
    gq:f32,mapq_avg:f32,baseq_avg:f32,
    vaf:f32,strand_bias:f32,
}

/// Caller stats
#[derive(Serialize, Deserialize, Debug, Default)] struct CallerStats {
    total_targets:usize, targets_with_variants:usize,
    total_variants_called:usize, variants_by_type:HashMap<String,usize>,
    elapsed_seconds:f64, threads_used:usize, params:HashMap<String,String>
}

fn main()->Result<()> {
    let start=Instant::now(); let cli=Cli::parse();
    let level = if cli.verbose { Level::DEBUG } else { Level::INFO };
    tracing_subscriber::fmt().with_max_level(level).init();
    let threads = if cli.threads==0 { num_cpus::get() } else { cli.threads };
    info!("Threads {}",threads);
    validate_inputs(&cli)?;
    info!("Using simplified BAM stub");
    let mut reader=SimpleBamReader::new(&cli.bam)?;
    let regions=get_regions(&mut reader,&cli)?;
    info!("Regions {}",regions.len());
    let stats=Arc::new(Mutex::new(CallerStats{total_targets:regions.len(),threads_used:threads,..Default::default()}));
    let pb=ProgressBar::new(regions.len() as u64);
    pb.set_style(ProgressStyle::default_bar().template("{bar:40.cyan/blue} {pos}/{len}").unwrap());
    let mut all_calls=Vec::new();
    for region in &regions {
        debug!("Region {}",region.name);
        let calls=generate_mock_calls(region,10);
        if !calls.is_empty(){ let mut s=stats.lock().unwrap(); s.targets_with_variants+=1; s.total_variants_called+=calls.len(); for c in &calls{ let t= if is_transition(c.ref_base,c.alt_base) {"transition"} else {"other"}; *s.variants_by_type.entry(t.into()).or_insert(0)+=1;} all_calls.extend(calls);}        
        pb.inc(1);
    }
    pb.finish_with_message("done");
    if all_calls.is_empty(){ warn!("No variants"); return Err(CallerError::NoVariants.into()); }
    export_variants(&all_calls,&cli.out)?;
    info!("Exported {}",all_calls.len());
    let mut s=stats.lock().unwrap(); s.elapsed_seconds=start.elapsed().as_secs_f64(); s.params.insert("min_depth".into(),cli.min_depth.to_string());
    if let Some(p)=&cli.stats{ export_stats(&s,p)?; info!("Stats at {}",p.display()); }
    print_summary(&all_calls,&s);
    Ok(())
}

fn validate_inputs(cli: &Cli) -> Result<()> {
    // Check if BAM file exists
    if !cli.bam.exists() {
        return Err(anyhow!("BAM file does not exist: {}", cli.bam.display()));
    }
    // Check FASTA if provided
    if let Some(fasta_path) = &cli.fasta {
        if !fasta_path.exists() {
            return Err(anyhow!("FASTA file does not exist: {}", fasta_path.display()));
        }
    }
    Ok(())
}

fn get_regions(reader: &mut SimpleBamReader, cli: &Cli) -> Result<Vec<Region>> {
    // Read header for reference sequences
    let header = reader.read_header()?;
    let refs = header.reference_sequences();

    // If a region is specified, parse and validate it
    if let Some(region_str) = &cli.region {
        let region = Region::from_str(region_str)?;
        if !refs.contains_key(&region.name) {
            return Err(anyhow!("Region '{}' not found in BAM header", region.name));
        }
        return Ok(vec![region]);
    }

    // Otherwise, return all contigs
    let regions = refs
        .keys()
        .map(|name| Region { name: name.clone(), start: None, end: None })
        .collect();
    Ok(regions)
}

fn generate_mock_calls(region: &Region, count: usize) -> Vec<Call> {
    let mut calls = Vec::with_capacity(count);
    for i in 0..count {
        let pos = region.start.unwrap_or(1000) as i64 + (i as i64) * 100;
        calls.push(Call {
            chrom: region.name.clone(),
            pos,
            ref_base: 'A',
            alt_base: 'C',
            depth: 30,
            ref_count: 20,
            alt_count: 10,
            gq: 30.0,
            mapq_avg: 40.0,
            baseq_avg: 35.0,
            vaf: 0.33,
            strand_bias: 0.9,
        });
    }
    calls
}

fn is_transition(r: char, a: char) -> bool {
    matches!(
        (r.to_ascii_uppercase(), a.to_ascii_uppercase()),
        ('A', 'G') | ('G', 'A') | ('C', 'T') | ('T', 'C')
    )
}

fn is_transversion(r: char, a: char) -> bool {
    let r = r.to_ascii_uppercase();
    let a = a.to_ascii_uppercase();
    (r == 'A' || r == 'G') && (a == 'C' || a == 'T') ||
    (r == 'C' || r == 'T') && (a == 'A' || a == 'G')
}

fn export_variants(calls: &[Call], out: &Path) -> Result<()> {
    // Convert calls to DataFrame
    let mut df = calls_to_dataframe(calls)?;

    // Write DataFrame to Parquet file
    let file = File::create(out).context("creating output file failed")?;
    let mut writer = BufWriter::new(file);
    ParquetWriter::new(&mut writer)
        .finish(&mut df)
        .context("writing Parquet file failed")?;
    Ok(())
}

fn calls_to_dataframe(calls: &[Call]) -> Result<DataFrame> {
    if calls.is_empty() {
        return Err(anyhow!("No variants to convert to DataFrame"));
    }
    let chrom = Series::new(
        "chrom",
        calls.iter().map(|c| c.chrom.clone()).collect::<Vec<String>>(),
    );
    let pos = Series::new(
        "pos",
        calls.iter().map(|c| c.pos).collect::<Vec<i64>>(),
    );
    DataFrame::new(vec![chrom, pos]).context("failed to create DataFrame")
}

fn export_stats(stats: &CallerStats, out: &Path) -> Result<()> {
    let json = serde_json::to_string_pretty(stats)?;
    std::fs::write(out, json).context("writing stats JSON failed")?;
    Ok(())
}

fn print_summary(calls: &[Call], stats: &CallerStats) {
    println!("=== Variant Calling Summary ===");
    println!("Total targets processed: {}", stats.total_targets);
    println!("Targets with variants: {}", stats.targets_with_variants);
    println!("Total variants called: {}", stats.total_variants_called);
    println!("Transition/Transversion ratio: {:.2}", 
        stats.variants_by_type.get("transition").cloned().unwrap_or(0) as f64
        / stats.variants_by_type.get("transversion").cloned().unwrap_or(1) as f64
    );
    println!("Runtime: {:.2} seconds", stats.elapsed_seconds);
    println!("Threads used: {}", stats.threads_used);
    println!("Parameters:");
    for (k, v) in &stats.params {
        println!("  {}: {}", k, v);
    }
}
