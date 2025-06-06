use std::collections::HashMap;
use std::path::PathBuf;
use std::fs::File;
use std::io::Write;
use anyhow::{Context, Result, bail};
use clap::Parser;
use log::{info, warn, error, debug};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use statrs::correlation::pearson;
use statrs::statistics::{Statistics, OrderStatistics};
use ndarray::{Array1, Array2, Axis};
use ndarray_stats::CorrelationExt;
use itertools::Itertools;
use indicatif::{ProgressBar, ProgressStyle};
use chrono::{DateTime, Utc};
use polars::prelude::*;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Coverage data file (CSV format)
    #[arg(short, long)]
    coverage: PathBuf,
    
    /// Expression data file (CSV format)
    #[arg(short, long)]
    expression: PathBuf,
    
    /// Output integration results file (JSON format)
    #[arg(short, long, default_value = "integration_results.json")]
    output: PathBuf,
    
    /// Gene annotation file (GTF/GFF format, optional)
    #[arg(long)]
    annotation: Option<PathBuf>,
    
    /// Sample metadata file (CSV format, optional)
    #[arg(long)]
    metadata: Option<PathBuf>,
    
    /// Correlation method (pearson, spearman, kendall)
    #[arg(long, default_value = "pearson")]
    correlation_method: String,
    
    /// Minimum correlation threshold for reporting
    #[arg(long, default_value_t = 0.3)]
    min_correlation: f64,
    
    /// P-value threshold for significance
    #[arg(long, default_value_t = 0.05)]
    pvalue_threshold: f64,
    
    /// Minimum number of samples for correlation
    #[arg(long, default_value_t = 3)]
    min_samples: usize,
    
    /// Output detailed correlation matrix
    #[arg(long)]
    output_matrix: Option<PathBuf>,
    
    /// Output network file (GraphML format)
    #[arg(long)]
    output_network: Option<PathBuf>,
    
    /// Generate visualization plots
    #[arg(long)]
    generate_plots: bool,
    
    /// Number of threads to use
    #[arg(short, long, default_value_t = 4)]
    threads: usize,
    
    /// Enable verbose logging
    #[arg(short, long)]
    verbose: bool,
}

#[derive(Debug, Serialize, Deserialize)]
struct IntegrationResult {
    gene_id: String,
    gene_name: Option<String>,
    gene_type: Option<String>,
    chromosome: Option<String>,
    start: Option<u64>,
    end: Option<u64>,
    coverage_stats: CoverageStats,
    expression_stats: ExpressionStats,
    correlation_analysis: CorrelationAnalysis,
    multi_omics_score: f64,
}

#[derive(Debug, Serialize, Deserialize)]
struct CoverageStats {
    mean_coverage: f64,
    median_coverage: f64,
    max_coverage: f64,
    coverage_variance: f64,
    samples_with_coverage: usize,
    peak_count: Option<usize>,
}

#[derive(Debug, Serialize, Deserialize)]
struct ExpressionStats {
    mean_expression: f64,
    median_expression: f64,
    max_expression: f64,
    expression_variance: f64,
    samples_with_expression: usize,
    fold_change: Option<f64>,
}

#[derive(Debug, Serialize, Deserialize)]
struct CorrelationAnalysis {
    correlation_coefficient: f64,
    pvalue: Option<f64>,
    confidence_interval: Option<(f64, f64)>,
    sample_count: usize,
    correlation_method: String,
    significance_level: String,
}

#[derive(Debug, Serialize, Deserialize)]
struct IntegrationSummary {
    analysis_date: DateTime<Utc>,
    total_genes: usize,
    significant_correlations: usize,
    positive_correlations: usize,
    negative_correlations: usize,
    mean_correlation: f64,
    correlation_distribution: CorrelationDistribution,
    quality_metrics: QualityMetrics,
    processing_stats: ProcessingStats,
}

#[derive(Debug, Serialize, Deserialize)]
struct CorrelationDistribution {
    strong_positive: usize,    // > 0.7
    moderate_positive: usize,  // 0.3 - 0.7
    weak_positive: usize,      // 0.1 - 0.3
    weak_negative: usize,      // -0.3 - -0.1
    moderate_negative: usize,  // -0.7 - -0.3
    strong_negative: usize,    // < -0.7
}

#[derive(Debug, Serialize, Deserialize)]
struct QualityMetrics {
    data_completeness: f64,
    coverage_dynamic_range: f64,
    expression_dynamic_range: f64,
    sample_correlation_consistency: f64,
}

#[derive(Debug, Serialize, Deserialize)]
struct ProcessingStats {
    analysis_time_seconds: f64,
    genes_processed: usize,
    samples_analyzed: usize,
    memory_usage_mb: f64,
    correlation_method_used: String,
}

#[derive(Debug)]
struct MultiOmicsData {
    coverage_data: DataFrame,
    expression_data: DataFrame,
    gene_annotations: Option<HashMap<String, GeneAnnotation>>,
    sample_metadata: Option<HashMap<String, SampleMetadata>>,
}

#[derive(Debug, Clone)]
struct GeneAnnotation {
    gene_name: String,
    gene_type: String,
    chromosome: String,
    start: u64,
    end: u64,
    strand: String,
}

#[derive(Debug, Clone)]
struct SampleMetadata {
    sample_type: String,
    condition: String,
    batch: String,
    additional_info: HashMap<String, String>,
}

struct MultiOmicsIntegrator {
    args: Args,
}

impl MultiOmicsIntegrator {
    fn new(args: Args) -> Self {
        Self { args }
    }
    
    fn load_data(&self) -> Result<MultiOmicsData> {
        info!("Loading multi-omics data...");
        
        // Load coverage data
        info!("Loading coverage data from {:?}", self.args.coverage);
        let coverage_data = LazyFrame::scan_csv(&self.args.coverage, ScanArgsCSV::default())
            .context("Failed to load coverage data")?
            .collect()
            .context("Failed to collect coverage data")?;
        
        info!("Loaded coverage data: {} genes, {} samples", 
              coverage_data.height(), coverage_data.width() - 1);
        
        // Load expression data
        info!("Loading expression data from {:?}", self.args.expression);
        let expression_data = LazyFrame::scan_csv(&self.args.expression, ScanArgsCSV::default())
            .context("Failed to load expression data")?
            .collect()
            .context("Failed to collect expression data")?;
        
        info!("Loaded expression data: {} genes, {} samples", 
              expression_data.height(), expression_data.width() - 1);
        
        // Load gene annotations if provided
        let gene_annotations = if let Some(ref annotation_path) = self.args.annotation {
            info!("Loading gene annotations from {:?}", annotation_path);
            Some(self.load_gene_annotations(annotation_path)?)
        } else {
            None
        };
        
        // Load sample metadata if provided
        let sample_metadata = if let Some(ref metadata_path) = self.args.metadata {
            info!("Loading sample metadata from {:?}", metadata_path);
            Some(self.load_sample_metadata(metadata_path)?)
        } else {
            None
        };
        
        Ok(MultiOmicsData {
            coverage_data,
            expression_data,
            gene_annotations,
            sample_metadata,
        })
    }
    
    fn load_gene_annotations(&self, annotation_path: &PathBuf) -> Result<HashMap<String, GeneAnnotation>> {
        let mut annotations = HashMap::new();
        
        // Simple GTF/GFF parser (simplified for demonstration)
        let content = std::fs::read_to_string(annotation_path)?;
        for line in content.lines() {
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 9 && fields[2] == "gene" {
                let chromosome = fields[0].to_string();
                let start: u64 = fields[3].parse().unwrap_or(0);
                let end: u64 = fields[4].parse().unwrap_or(0);
                let strand = fields[6].to_string();
                
                // Parse attributes (simplified)
                let attributes = fields[8];
                let gene_id = self.extract_attribute(attributes, "gene_id").unwrap_or_default();
                let gene_name = self.extract_attribute(attributes, "gene_name").unwrap_or(gene_id.clone());
                let gene_type = self.extract_attribute(attributes, "gene_type").unwrap_or("unknown".to_string());
                
                annotations.insert(gene_id.clone(), GeneAnnotation {
                    gene_name,
                    gene_type,
                    chromosome,
                    start,
                    end,
                    strand,
                });
            }
        }
        
        info!("Loaded {} gene annotations", annotations.len());
        Ok(annotations)
    }
    
    fn extract_attribute(&self, attributes: &str, key: &str) -> Option<String> {
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
    
    fn load_sample_metadata(&self, metadata_path: &PathBuf) -> Result<HashMap<String, SampleMetadata>> {
        let df = LazyFrame::scan_csv(metadata_path, ScanArgsCSV::default())
            .context("Failed to load sample metadata")?
            .collect()?;
        
        let mut metadata = HashMap::new();
        
        // Assuming first column is sample_id
        let sample_ids = df.column("sample_id")
            .or_else(|_| df.get_column_by_index(0))
            .context("No sample_id column found")?;
        
        for i in 0..df.height() {
            if let Ok(sample_id) = sample_ids.get(i).unwrap().try_extract::<String>() {
                let sample_type = df.column("sample_type")
                    .and_then(|col| col.get(i))
                    .and_then(|val| val.try_extract::<String>().ok())
                    .unwrap_or_default();
                
                let condition = df.column("condition")
                    .and_then(|col| col.get(i))
                    .and_then(|val| val.try_extract::<String>().ok())
                    .unwrap_or_default();
                
                let batch = df.column("batch")
                    .and_then(|col| col.get(i))
                    .and_then(|val| val.try_extract::<String>().ok())
                    .unwrap_or_default();
                
                metadata.insert(sample_id, SampleMetadata {
                    sample_type,
                    condition,
                    batch,
                    additional_info: HashMap::new(),
                });
            }
        }
        
        info!("Loaded metadata for {} samples", metadata.len());
        Ok(metadata)
    }
    
    fn integrate_omics_data(&self, data: &MultiOmicsData) -> Result<Vec<IntegrationResult>> {
        info!("Starting multi-omics integration analysis...");
        
        // Get common genes between coverage and expression data
        let coverage_genes = self.get_gene_ids(&data.coverage_data)?;
        let expression_genes = self.get_gene_ids(&data.expression_data)?;
        
        let common_genes: Vec<String> = coverage_genes
            .intersection(&expression_genes)
            .cloned()
            .collect();
        
        info!("Found {} common genes between datasets", common_genes.len());
        
        if common_genes.is_empty() {
            bail!("No common genes found between coverage and expression datasets");
        }
        
        // Set up progress tracking
        let pb = ProgressBar::new(common_genes.len() as u64);
        pb.set_style(ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos:>7}/{len:7} {msg}")
            .unwrap());
        
        // Process genes in parallel
        let results: Result<Vec<_>> = common_genes
            .par_iter()
            .map(|gene_id| {
                pb.inc(1);
                self.analyze_gene(gene_id, data)
            })
            .collect();
        
        pb.finish_with_message("Integration analysis completed");
        
        let mut integration_results = results?;
        
        // Filter by correlation threshold
        integration_results.retain(|result| {
            result.correlation_analysis.correlation_coefficient.abs() >= self.args.min_correlation
        });
        
        // Sort by correlation strength
        integration_results.sort_by(|a, b| {
            b.correlation_analysis.correlation_coefficient
                .abs()
                .partial_cmp(&a.correlation_analysis.correlation_coefficient.abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        
        info!("Found {} genes with significant correlations", integration_results.len());
        Ok(integration_results)
    }
    
    fn get_gene_ids(&self, df: &DataFrame) -> Result<std::collections::HashSet<String>> {
        let gene_column = df.column("gene_id")
            .or_else(|_| df.get_column_by_index(0))
            .context("No gene_id column found")?;
        
        let mut genes = std::collections::HashSet::new();
        for i in 0..df.height() {
            if let Ok(gene_id) = gene_column.get(i).unwrap().try_extract::<String>() {
                genes.insert(gene_id);
            }
        }
        
        Ok(genes)
    }
    
    fn analyze_gene(&self, gene_id: &str, data: &MultiOmicsData) -> Result<IntegrationResult> {
        // Extract coverage values for this gene
        let coverage_values = self.extract_gene_values(gene_id, &data.coverage_data)?;
        let expression_values = self.extract_gene_values(gene_id, &data.expression_data)?;
        
        // Ensure we have enough samples
        if coverage_values.len() < self.args.min_samples || expression_values.len() < self.args.min_samples {
            bail!("Insufficient samples for gene {}", gene_id);
        }
        
        // Calculate statistics
        let coverage_stats = self.calculate_coverage_stats(&coverage_values);
        let expression_stats = self.calculate_expression_stats(&expression_values);
        
        // Perform correlation analysis
        let correlation_analysis = self.perform_correlation_analysis(&coverage_values, &expression_values)?;
        
        // Calculate multi-omics integration score
        let multi_omics_score = self.calculate_integration_score(&correlation_analysis, &coverage_stats, &expression_stats);
        
        // Get gene annotation if available
        let annotation = data.gene_annotations.as_ref().and_then(|ann| ann.get(gene_id));
        
        Ok(IntegrationResult {
            gene_id: gene_id.to_string(),
            gene_name: annotation.map(|a| a.gene_name.clone()),
            gene_type: annotation.map(|a| a.gene_type.clone()),
            chromosome: annotation.map(|a| a.chromosome.clone()),
            start: annotation.map(|a| a.start),
            end: annotation.map(|a| a.end),
            coverage_stats,
            expression_stats,
            correlation_analysis,
            multi_omics_score,
        })
    }
    
    fn extract_gene_values(&self, gene_id: &str, df: &DataFrame) -> Result<Vec<f64>> {
        // Find the row for this gene
        let gene_column = df.column("gene_id")
            .or_else(|_| df.get_column_by_index(0))
            .context("No gene_id column found")?;
        
        let mut gene_row_idx = None;
        for i in 0..df.height() {
            if let Ok(current_gene) = gene_column.get(i).unwrap().try_extract::<String>() {
                if current_gene == gene_id {
                    gene_row_idx = Some(i);
                    break;
                }
            }
        }
        
        let row_idx = gene_row_idx.context(format!("Gene {} not found", gene_id))?;
        
        // Extract values from all sample columns (skip first column which is gene_id)
        let mut values = Vec::new();
        for col_idx in 1..df.width() {
            if let Ok(column) = df.get_column_by_index(col_idx) {
                if let Ok(value) = column.get(row_idx).unwrap().try_extract::<f64>() {
                    if !value.is_nan() && value.is_finite() {
                        values.push(value);
                    }
                }
            }
        }
        
        Ok(values)
    }
    
    fn calculate_coverage_stats(&self, values: &[f64]) -> CoverageStats {
        if values.is_empty() {
            return CoverageStats {
                mean_coverage: 0.0,
                median_coverage: 0.0,
                max_coverage: 0.0,
                coverage_variance: 0.0,
                samples_with_coverage: 0,
                peak_count: None,
            };
        }
        
        let mean_coverage = values.iter().sum::<f64>() / values.len() as f64;
        let mut sorted_values = values.to_vec();
        sorted_values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let median_coverage = sorted_values[sorted_values.len() / 2];
        let max_coverage = sorted_values.last().copied().unwrap_or(0.0);
        
        let variance = values.iter()
            .map(|x| (x - mean_coverage).powi(2))
            .sum::<f64>() / values.len() as f64;
        
        let samples_with_coverage = values.iter().filter(|&&x| x > 0.0).count();
        
        CoverageStats {
            mean_coverage,
            median_coverage,
            max_coverage,
            coverage_variance: variance,
            samples_with_coverage,
            peak_count: None, // Would be calculated from peak data if available
        }
    }
    
    fn calculate_expression_stats(&self, values: &[f64]) -> ExpressionStats {
        if values.is_empty() {
            return ExpressionStats {
                mean_expression: 0.0,
                median_expression: 0.0,
                max_expression: 0.0,
                expression_variance: 0.0,
                samples_with_expression: 0,
                fold_change: None,
            };
        }
        
        let mean_expression = values.iter().sum::<f64>() / values.len() as f64;
        let mut sorted_values = values.to_vec();
        sorted_values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let median_expression = sorted_values[sorted_values.len() / 2];
        let max_expression = sorted_values.last().copied().unwrap_or(0.0);
        
        let variance = values.iter()
            .map(|x| (x - mean_expression).powi(2))
            .sum::<f64>() / values.len() as f64;
        
        let samples_with_expression = values.iter().filter(|&&x| x > 0.0).count();
        
        ExpressionStats {
            mean_expression,
            median_expression,
            max_expression,
            expression_variance: variance,
            samples_with_expression,
            fold_change: None, // Would be calculated with control/treatment comparison
        }
    }
    
    fn perform_correlation_analysis(&self, coverage: &[f64], expression: &[f64]) -> Result<CorrelationAnalysis> {
        if coverage.len() != expression.len() {
            bail!("Coverage and expression vectors must have the same length");
        }
        
        if coverage.len() < 3 {
            bail!("Need at least 3 samples for correlation analysis");
        }
        
        // Calculate correlation coefficient based on method
        let correlation_coefficient = match self.args.correlation_method.as_str() {
            "pearson" => self.calculate_pearson_correlation(coverage, expression)?,
            "spearman" => self.calculate_spearman_correlation(coverage, expression)?,
            "kendall" => self.calculate_kendall_correlation(coverage, expression)?,
            _ => bail!("Unsupported correlation method: {}", self.args.correlation_method),
        };
        
        // Calculate p-value (simplified implementation)
        let pvalue = self.calculate_correlation_pvalue(correlation_coefficient, coverage.len());
        
        // Determine significance level
        let significance_level = if let Some(p) = pvalue {
            if p < 0.001 {
                "***".to_string()
            } else if p < 0.01 {
                "**".to_string()
            } else if p < 0.05 {
                "*".to_string()
            } else {
                "ns".to_string()
            }
        } else {
            "unknown".to_string()
        };
        
        // Calculate confidence interval (simplified)
        let confidence_interval = self.calculate_confidence_interval(correlation_coefficient, coverage.len());
        
        Ok(CorrelationAnalysis {
            correlation_coefficient,
            pvalue,
            confidence_interval,
            sample_count: coverage.len(),
            correlation_method: self.args.correlation_method.clone(),
            significance_level,
        })
    }
    
    fn calculate_pearson_correlation(&self, x: &[f64], y: &[f64]) -> Result<f64> {
        if x.len() != y.len() || x.is_empty() {
            bail!("Invalid input for correlation calculation");
        }
        
        let n = x.len() as f64;
        let sum_x: f64 = x.iter().sum();
        let sum_y: f64 = y.iter().sum();
        let sum_xx: f64 = x.iter().map(|val| val * val).sum();
        let sum_yy: f64 = y.iter().map(|val| val * val).sum();
        let sum_xy: f64 = x.iter().zip(y.iter()).map(|(a, b)| a * b).sum();
        
        let numerator = n * sum_xy - sum_x * sum_y;
        let denominator = ((n * sum_xx - sum_x * sum_x) * (n * sum_yy - sum_y * sum_y)).sqrt();
        
        if denominator == 0.0 {
            Ok(0.0)
        } else {
            Ok(numerator / denominator)
        }
    }
    
    fn calculate_spearman_correlation(&self, x: &[f64], y: &[f64]) -> Result<f64> {
        // Convert to ranks and calculate Pearson correlation of ranks
        let x_ranks = self.rank_transform(x);
        let y_ranks = self.rank_transform(y);
        self.calculate_pearson_correlation(&x_ranks, &y_ranks)
    }
    
    fn calculate_kendall_correlation(&self, x: &[f64], y: &[f64]) -> Result<f64> {
        // Simplified Kendall's tau implementation
        let n = x.len();
        let mut concordant = 0;
        let mut discordant = 0;
        
        for i in 0..n {
            for j in (i + 1)..n {
                let x_diff = x[i] - x[j];
                let y_diff = y[i] - y[j];
                
                if (x_diff > 0.0 && y_diff > 0.0) || (x_diff < 0.0 && y_diff < 0.0) {
                    concordant += 1;
                } else if (x_diff > 0.0 && y_diff < 0.0) || (x_diff < 0.0 && y_diff > 0.0) {
                    discordant += 1;
                }
            }
        }
        
        let total_pairs = n * (n - 1) / 2;
        if total_pairs == 0 {
            Ok(0.0)
        } else {
            Ok((concordant as f64 - discordant as f64) / total_pairs as f64)
        }
    }
    
    fn rank_transform(&self, values: &[f64]) -> Vec<f64> {
        let mut indexed_values: Vec<(usize, f64)> = values.iter().copied().enumerate().collect();
        indexed_values.sort_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        
        let mut ranks = vec![0.0; values.len()];
        for (rank, (original_index, _)) in indexed_values.iter().enumerate() {
            ranks[*original_index] = (rank + 1) as f64;
        }
        
        ranks
    }
    
    fn calculate_correlation_pvalue(&self, r: f64, n: usize) -> Option<f64> {
        if n < 3 {
            return None;
        }
        
        // Simplified t-test for correlation significance
        let df = n as f64 - 2.0;
        let t_stat = r * (df / (1.0 - r * r)).sqrt();
        
        // Very simplified p-value calculation (in practice, use proper statistical library)
        let p_value = if t_stat.abs() > 2.576 {
            0.01
        } else if t_stat.abs() > 1.96 {
            0.05
        } else {
            0.1
        };
        
        Some(p_value)
    }
    
    fn calculate_confidence_interval(&self, r: f64, n: usize) -> Option<(f64, f64)> {
        if n < 4 {
            return None;
        }
        
        // Fisher's z-transformation for confidence interval
        let z = 0.5 * ((1.0 + r) / (1.0 - r)).ln();
        let