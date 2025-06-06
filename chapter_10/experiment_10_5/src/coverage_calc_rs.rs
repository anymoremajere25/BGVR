use std::collections::HashMap;
use std::path::PathBuf;
use std::fs::File;
use std::io::{Write, BufReader, BufRead};
use anyhow::{Context, Result, bail};
use clap::Parser;
use log::{info, warn};
use rayon::prelude::*;
use rust_htslib::bam::{self, Read, HeaderView};
use serde::{Deserialize, Serialize};

#[derive(Parser, Debug)]
#[command(author, version, about = "Calculate coverage from BAM files for multi-omics integration")]
struct Args {
    /// Input BAM file
    bam_file: PathBuf,
    
    /// Gene annotation file (BED or GTF format)
    #[arg(short, long)]
    annotation: PathBuf,
    
    /// Output coverage file (CSV format)
    #[arg(short, long, default_value = "coverage.csv")]
    output: PathBuf,
    
    /// Coverage calculation method (total, mean, median, max)
    #[arg(short, long, default_value = "mean")]
    method: String,
    
    /// Extend gene regions by N bp (for promoter analysis)
    #[arg(long, default_value_t = 0)]
    extend_regions: u64,
    
    /// Minimum mapping quality
    #[arg(long, default_value_t = 10)]
    min_mapq: u8,
    
    /// Count only unique reads (remove duplicates)
    #[arg(long)]
    unique_only: bool,
    
    /// Normalize by gene length (RPKM-style)
    #[arg(long)]
    normalize_length: bool,
    
    /// Number of threads
    #[arg(short, long, default_value_t = 4)]
    threads: usize,
    
    /// Enable verbose logging
    #[arg(short, long)]
    verbose: bool,
}

#[derive(Debug, Clone)]
struct GeneRegion {
    gene_id: String,
    chromosome: String,
    start: u64,
    end: u64,
    strand: char,
    gene_name: Option<String>,
}

#[derive(Debug, Serialize)]
struct CoverageResult {
    gene_id: String,
    gene_name: String,
    chromosome: String,
    start: u64,
    end: u64,
    length: u64,
    raw_coverage: f64,
    normalized_coverage: f64,
    read_count: u64,
}

struct CoverageCalculator {
    args: Args,
    gene_regions: Vec<GeneRegion>,
}

impl CoverageCalculator {
    fn new(args: Args) -> Result<Self> {
        let gene_regions = Self::load_gene_regions(&args.annotation, args.extend_regions)?;
        info!("Loaded {} gene regions", gene_regions.len());
        
        Ok(Self {
            args,
            gene_regions,
        })
    }
    
    fn load_gene_regions(annotation_path: &PathBuf, extend_bp: u64) -> Result<Vec<GeneRegion>> {
        let file = File::open(annotation_path)
            .with_context(|| format!("Failed to open annotation file: {:?}", annotation_path))?;
        
        let reader = BufReader::new(file);
        let mut regions = Vec::new();
        
        // Detect file format
        let first_line = std::fs::read_to_string(annotation_path)?
            .lines()
            .find(|line| !line.starts_with('#') && !line.trim().is_empty())
            .unwrap_or("")
            .to_string();
        
        if first_line.split('\t').count() >= 9 {
            // GTF/GFF format
            regions = Self::parse_gtf_file(annotation_path, extend_bp)?;
        } else if first_line.split('\t').count() >= 6 {
            // BED format
            regions = Self::parse_bed_file(annotation_path, extend_bp)?;
        } else {
            bail!("Unsupported annotation file format. Expected GTF or BED format.");
        }
        
        Ok(regions)
    }
    
    fn parse_gtf_file(gtf_path: &PathBuf, extend_bp: u64) -> Result<Vec<GeneRegion>> {
        let content = std::fs::read_to_string(gtf_path)?;
        let mut regions = Vec::new();
        
        for line in content.lines() {
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 9 && fields[2] == "gene" {
                let chromosome = fields[0].to_string();
                let start: u64 = fields[3].parse().unwrap_or(0).saturating_sub(extend_bp);
                let end: u64 = fields[4].parse().unwrap_or(0) + extend_bp;
                let strand = fields[6].chars().next().unwrap_or('.');
                
                // Parse attributes
                let attributes = fields[8];
                let gene_id = Self::extract_gtf_attribute(attributes, "gene_id")
                    .unwrap_or_else(|| format!("unknown_{}", regions.len()));
                let gene_name = Self::extract_gtf_attribute(attributes, "gene_name");
                
                regions.push(GeneRegion {
                    gene_id,
                    chromosome,
                    start,
                    end,
                    strand,
                    gene_name,
                });
            }
        }
        
        Ok(regions)
    }
    
    fn parse_bed_file(bed_path: &PathBuf, extend_bp: u64) -> Result<Vec<GeneRegion>> {
        let content = std::fs::read_to_string(bed_path)?;
        let mut regions = Vec::new();
        
        for (line_num, line) in content.lines().enumerate() {
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 6 {
                let chromosome = fields[0].to_string();
                let start: u64 = fields[1].parse()
                    .with_context(|| format!("Invalid start coordinate on line {}", line_num + 1))?
                    .saturating_sub(extend_bp);
                let end: u64 = fields[2].parse()
                    .with_context(|| format!("Invalid end coordinate on line {}", line_num + 1))?
                    + extend_bp;
                
                let gene_id = if fields.len() > 3 && !fields[3].is_empty() {
                    fields[3].to_string()
                } else {
                    format!("gene_{}", line_num + 1)
                };
                
                let strand = if fields.len() > 5 {
                    fields[5].chars().next().unwrap_or('.')
                } else {
                    '.'
                };
                
                let gene_name = if fields.len() > 6 && !fields[6].is_empty() {
                    Some(fields[6].to_string())
                } else {
                    None
                };
                
                regions.push(GeneRegion {
                    gene_id,
                    chromosome,
                    start,
                    end,
                    strand,
                    gene_name,
                });
            }
        }
        
        Ok(regions)
    }
    
    fn extract_gtf_attribute(attributes: &str, key: &str) -> Option<String> {
        for attr in attributes.split(';') {
            let attr = attr.trim();
            if attr.starts_with(key) {
                if let Some(value) = attr.split('=').nth(1).or_else(|| attr.split(' ').nth(1)) {
                    return Some(value.trim_matches('"').to_string());
                }
            }
        }
        None
    }
    
    fn calculate_coverage(&self) -> Result<Vec<CoverageResult>> {
        info!("Calculating coverage for {} genes using method: {}", 
              self.gene_regions.len(), self.args.method);
        
        let mut bam_reader = bam::Reader::from_path(&self.args.bam_file)?;
        let header = bam_reader.header().clone();
        
        // Calculate coverage for each gene region in parallel
        let results: Result<Vec<_>> = self.gene_regions
            .par_iter()
            .map(|region| self.calculate_gene_coverage(region, &header))
            .collect();
        
        let mut coverage_results = results?;
        
        // Apply length normalization if requested
        if self.args.normalize_length {
            let total_reads = coverage_results.iter().map(|r| r.read_count).sum::<u64>() as f64;
            let total_length = coverage_results.iter().map(|r| r.length).sum::<u64>() as f64;
            
            for result in &mut coverage_results {
                // RPKM-style normalization: (reads * 1e9) / (gene_length * total_reads)
                if total_reads > 0.0 && result.length > 0 {
                    result.normalized_coverage = (result.read_count as f64 * 1e9) / 
                        (result.length as f64 * total_reads);
                } else {
                    result.normalized_coverage = 0.0;
                }
            }
        }
        
        // Sort by coverage (descending)
        coverage_results.sort_by(|a, b| 
            b.raw_coverage.partial_cmp(&a.raw_coverage).unwrap_or(std::cmp::Ordering::Equal));
        
        Ok(coverage_results)
    }
    
    fn calculate_gene_coverage(&self, region: &GeneRegion, header: &HeaderView) -> Result<CoverageResult> {
        // Open a new BAM reader for this thread
        let mut bam_reader = bam::Reader::from_path(&self.args.bam_file)?;
        
        // Get chromosome ID
        let tid = header.tid(region.chromosome.as_bytes());
        if tid.is_none() {
            warn!("Chromosome {} not found in BAM header", region.chromosome);
            return Ok(CoverageResult {
                gene_id: region.gene_id.clone(),
                gene_name: region.gene_name.clone().unwrap_or_else(|| region.gene_id.clone()),
                chromosome: region.chromosome.clone(),
                start: region.start,
                end: region.end,
                length: region.end - region.start,
                raw_coverage: 0.0,
                normalized_coverage: 0.0,
                read_count: 0,
            });
        }
        
        let tid = tid.unwrap();
        
        // Fetch reads in the region
        if let Err(e) = bam_reader.fetch((tid, region.start as i64, region.end as i64)) {
            warn!("Failed to fetch region {}:{}-{}: {}", 
                  region.chromosome, region.start, region.end, e);
            return Ok(CoverageResult {
                gene_id: region.gene_id.clone(),
                gene_name: region.gene_name.clone().unwrap_or_else(|| region.gene_id.clone()),
                chromosome: region.chromosome.clone(),
                start: region.start,
                end: region.end,
                length: region.end - region.start,
                raw_coverage: 0.0,
                normalized_coverage: 0.0,
                read_count: 0,
            });
        }
        
        // Collect coverage values
        let mut coverage_values = Vec::new();
        let mut read_count = 0u64;
        let mut processed_positions = std::collections::HashSet::new();
        
        for result in bam_reader.records() {
            let record = result?;
            
            // Apply filters
            if record.is_unmapped() || 
               record.mapq() < self.args.min_mapq ||
               (self.args.unique_only && record.is_duplicate()) {
                continue;
            }
            
            // Check if read overlaps with gene region
            let read_start = record.pos() as u64;
            let read_end = read_start + record.seq_len() as u64;
            
            if read_end >= region.start && read_start <= region.end {
                read_count += 1;
                
                // For position-based coverage, track each covered position
                for pos in read_start.max(region.start)..read_end.min(region.end) {
                    processed_positions.insert(pos);
                }
            }
        }
        
        // Calculate coverage based on method
        let raw_coverage = match self.args.method.as_str() {
            "total" => read_count as f64,
            "mean" => {
                let region_length = region.end - region.start;
                if region_length > 0 {
                    read_count as f64 / region_length as f64
                } else {
                    0.0
                }
            },
            "density" => {
                let covered_positions = processed_positions.len() as f64;
                let region_length = region.end - region.start;
                if region_length > 0 {
                    covered_positions / region_length as f64
                } else {
                    0.0
                }
            },
            _ => {
                // Default to mean
                let region_length = region.end - region.start;
                if region_length > 0 {
                    read_count as f64 / region_length as f64
                } else {
                    0.0
                }
            }
        };
        
        Ok(CoverageResult {
            gene_id: region.gene_id.clone(),
            gene_name: region.gene_name.clone().unwrap_or_else(|| region.gene_id.clone()),
            chromosome: region.chromosome.clone(),
            start: region.start,
            end: region.end,
            length: region.end - region.start,
            raw_coverage,
            normalized_coverage: raw_coverage, // Will be updated later if normalization is requested
            read_count,
        })
    }
    
    fn write_results(&self, results: &[CoverageResult]) -> Result<()> {
        info!("Writing coverage results to {:?}", self.args.output);
        
        let mut csv_content = String::new();
        csv_content.push_str("gene_id,gene_name,chromosome,start,end,length,raw_coverage,normalized_coverage,read_count\n");
        
        for result in results {
            csv_content.push_str(&format!(
                "{},{},{},{},{},{},{:.6},{:.6},{}\n",
                result.gene_id,
                result.gene_name,
                result.chromosome,
                result.start,
                result.end,
                result.length,
                result.raw_coverage,
                result.normalized_coverage,
                result.read_count
            ));
        }
        
        std::fs::write(&self.args.output, csv_content)
            .with_context(|| format!("Failed to write results to {:?}", self.args.output))?;
        
        Ok(())
    }
    
    pub fn run(&mut self) -> Result<()> {
        info!("Starting coverage calculation...");
        
        // Validate input files
        if !self.args.bam_file.exists() {
            bail!("BAM file does not exist: {:?}", self.args.bam_file);
        }
        
        if !self.args.annotation.exists() {
            bail!("Annotation file does not exist: {:?}", self.args.annotation);
        }
        
        // Calculate coverage
        let results = self.calculate_coverage()?;
        
        // Write results
        self.write_results(&results)?;
        
        // Print summary
        let total_genes = results.len();
        let genes_with_coverage = results.iter().filter(|r| r.raw_coverage > 0.0).count();
        let mean_coverage = if total_genes > 0 {
            results.iter().map(|r| r.raw_coverage).sum::<f64>() / total_genes as f64
        } else {
            0.0
        };
        let max_coverage = results.iter().map(|r| r.raw_coverage).fold(0.0, f64::max);
        
        println!("\n=== Coverage Calculation Summary ===");
        println!("Total genes: {}", total_genes);
        println!("Genes with coverage: {} ({:.1}%)", 
                 genes_with_coverage, 
                 100.0 * genes_with_coverage as f64 / total_genes as f64);
        println!("Mean coverage: {:.3}", mean_coverage);
        println!("Max coverage: {:.3}", max_coverage);
        println!("Method: {}", self.args.method);
        println!("Length normalized: {}", self.args.normalize_length);
        println!("===================================");
        
        info!("Coverage calculation completed successfully");
        Ok(())
    }
}

fn main() -> Result<()> {
    let args = Args::parse();
    
    // Initialize logging
    let log_level = if args.verbose { "debug" } else { "info" };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(log_level)).init();
    
    info!("Coverage calculator starting...");
    
    // Set thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .context("Failed to initialize thread pool")?;
    
    // Run coverage calculation
    let mut calculator = CoverageCalculator::new(args)?;
    calculator.run()?;
    
    Ok(())
}