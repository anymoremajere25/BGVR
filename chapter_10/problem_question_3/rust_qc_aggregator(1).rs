// src/bin/qc_aggregator.rs
use clap::Parser;
use serde::{Deserialize, Serialize};
use serde_json::json;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

#[derive(Parser)]
#[command(name = "qc_aggregator")]
#[command(about = "Aggregate ATAC-Seq QC metrics into comprehensive report")]
pub struct Args {
    #[arg(long, help = "Comma-separated list of TSS enrichment JSON files")]
    tss_files: String,
    
    #[arg(long, help = "File containing list of alignment stats JSON files")]
    alignment_stats: String,
    
    #[arg(long, help = "Output HTML report file")]
    output_html: String,
    
    #[arg(long, help = "Output JSON summary file")]
    output_json: String,
    
    #[arg(long, help = "Output TSV metrics table")]
    output_table: String,
    
    #[arg(long, help = "Pipeline version")]
    pipeline_version: String,
    
    #[arg(long, help = "Analysis date")]
    analysis_date: String,
    
    #[arg(long, help = "Log file path")]
    log_file: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct AlignmentStats {
    sample: String,
    total_reads: u64,
    mapped_reads: u64,
    high_quality_reads: u64,
    duplicate_reads: u64,
    mitochondrial_reads: u64,
    final_reads: u64,
    mapping_rate: f64,
    quality_rate: f64,
    duplicate_rate: f64,
    alignment_time_seconds: f64,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct TSSStats {
    sample: String,
    tss_enrichment: f64,
    signal_to_noise: f64,
    max_signal: f64,
    background_signal: f64,
    total_tss_sites: usize,
    window_size: i32,
    quality_grade: String,
    analysis_timestamp: String,
}

#[derive(Debug, Serialize)]
pub struct QCSummary {
    pipeline_info: PipelineInfo,
    sample_count: usize,
    overall_metrics: OverallMetrics,
    sample_metrics: Vec<SampleMetrics>,
    quality_thresholds: QualityThresholds,
    recommendations: Vec<String>,
}

#[derive(Debug, Serialize)]
pub struct PipelineInfo {
    version: String,
    analysis_date: String,
    total_samples: usize,
    successful_samples: usize,
}

#[derive(Debug, Serialize)]
pub struct OverallMetrics {
    avg_mapping_rate: f64,
    avg_tss_enrichment: f64,
    avg_signal_to_noise: f64,
    avg_duplicate_rate: f64,
    samples_passed_qc: usize,
}

#[derive(Debug, Serialize)]
pub struct SampleMetrics {
    sample_name: String,
    mapping_rate: f64,
    tss_enrichment: f64,
    signal_to_noise: f64,
    duplicate_rate: f64,
    final_reads: u64,
    quality_grade: String,
    qc_pass: bool,
}

#[derive(Debug, Serialize)]
pub struct QualityThresholds {
    min_mapping_rate: f64,
    min_tss_enrichment: f64,
    min_signal_to_noise: f64,
    max_duplicate_rate: f64,
    min_final_reads: u64,
}

impl QualityThresholds {
    fn default() -> Self {
        Self {
            min_mapping_rate: 70.0,
            min_tss_enrichment: 5.0,
            min_signal_to_noise: 3.0,
            max_duplicate_rate: 30.0,
            min_final_reads: 10_000_000,
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
    
    log_message(&mut log_writer, "Starting QC aggregation")?;
    log_message(&mut log_writer, &format!("Pipeline version: {}", args.pipeline_version))?;
    log_message(&mut log_writer, &format!("Analysis date: {}", args.analysis_date))?;
    
    // Step 1: Load TSS enrichment data
    log_message(&mut log_writer, "Loading TSS enrichment data...")?;
    let tss_files: Vec<&str> = args.tss_files.split(',').collect();
    let mut tss_data = HashMap::new();
    
    for tss_file in tss_files {
        if let Ok(data) = load_tss_data(tss_file.trim()) {
            tss_data.insert(data.sample.clone(), data);
            log_message(&mut log_writer, &format!("Loaded TSS data for sample: {}", data.sample))?;
        } else {
            log_message(&mut log_writer, &format!("Warning: Could not load TSS data from: {}", tss_file))?;
        }
    }
    
    // Step 2: Load alignment statistics
    log_message(&mut log_writer, "Loading alignment statistics...")?;
    let alignment_files = load_alignment_file_list(&args.alignment_stats)?;
    let mut alignment_data = HashMap::new();
    
    for alignment_file in alignment_files {
        if let Ok(data) = load_alignment_data(&alignment_file) {
            alignment_data.insert(data.sample.clone(), data);
            log_message(&mut log_writer, &format!("Loaded alignment data for sample: {}", data.sample))?;
        } else {
            log_message(&mut log_writer, &format!("Warning: Could not load alignment data from: {}", alignment_file))?;
        }
    }
    
    // Step 3: Merge data and create QC summary
    log_message(&mut log_writer, "Creating QC summary...")?;
    let qc_summary = create_qc_summary(
        &tss_data,
        &alignment_data,
        &args.pipeline_version,
        &args.analysis_date,
    )?;
    
    log_message(&mut log_writer, &format!("Total samples processed: {}", qc_summary.sample_count))?;
    log_message(&mut log_writer, &format!("Samples passed QC: {}", qc_summary.overall_metrics.samples_passed_qc))?;
    log_message(&mut log_writer, &format!("Average TSS enrichment: {:.2}", qc_summary.overall_metrics.avg_tss_enrichment))?;
    
    // Step 4: Generate outputs
    log_message(&mut log_writer, "Generating output files...")?;
    
    // Write JSON summary
    let json_output = serde_json::to_string_pretty(&qc_summary)?;
    std::fs::write(&args.output_json, json_output)?;
    log_message(&mut log_writer, &format!("JSON summary written to: {}", args.output_json))?;
    
    // Write TSV table
    write_metrics_table(&qc_summary, &args.output_table, &mut log_writer)?;
    
    // Write HTML report
    write_html_report(&qc_summary, &args.output_html, &mut log_writer)?;
    
    log_message(&mut log_writer, "QC aggregation completed successfully")?;
    
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

fn load_tss_data(file_path: &str) -> Result<TSSStats, Box<dyn std::error::Error>> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let tss_stats: TSSStats = serde_json::from_reader(reader)?;
    Ok(tss_stats)
}

fn load_alignment_data(file_path: &str) -> Result<AlignmentStats, Box<dyn std::error::Error>> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let alignment_stats: AlignmentStats = serde_json::from_reader(reader)?;
    Ok(alignment_stats)
}

fn load_alignment_file_list(list_file: &str) -> Result<Vec<String>, Box<dyn std::error::Error>> {
    let file = File::open(list_file)?;
    let reader = BufReader::new(file);
    let mut files = Vec::new();
    
    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if !trimmed.is_empty() && !trimmed.starts_with('#') {
            files.push(trimmed.to_string());
        }
    }
    
    Ok(files)
}

fn create_qc_summary(
    tss_data: &HashMap<String, TSSStats>,
    alignment_data: &HashMap<String, AlignmentStats>,
    pipeline_version: &str,
    analysis_date: &str,
) -> Result<QCSummary, Box<dyn std::error::Error>> {
    let thresholds = QualityThresholds::default();
    let mut sample_metrics = Vec::new();
    let mut total_mapping_rate = 0.0;
    let mut total_tss_enrichment = 0.0;
    let mut total_signal_to_noise = 0.0;
    let mut total_duplicate_rate = 0.0;
    let mut passed_qc_count = 0;
    
    // Get all unique sample names
    let mut all_samples = std::collections::HashSet::new();
    for sample in tss_data.keys() {
        all_samples.insert(sample.clone());
    }
    for sample in alignment_data.keys() {
        all_samples.insert(sample.clone());
    }
    
    for sample_name in all_samples {
        let tss = tss_data.get(&sample_name);
        let alignment = alignment_data.get(&sample_name);
        
        let mapping_rate = alignment.map(|a| a.mapping_rate).unwrap_or(0.0);
        let tss_enrichment = tss.map(|t| t.tss_enrichment).unwrap_or(0.0);
        let signal_to_noise = tss.map(|t| t.signal_to_noise).unwrap_or(0.0);
        let duplicate_rate = alignment.map(|a| a.duplicate_rate).unwrap_or(0.0);
        let final_reads = alignment.map(|a| a.final_reads).unwrap_or(0);
        let quality_grade = tss.map(|t| t.quality_grade.clone()).unwrap_or_else(|| "Unknown".to_string());
        
        // Determine if sample passes QC
        let qc_pass = mapping_rate >= thresholds.min_mapping_rate
            && tss_enrichment >= thresholds.min_tss_enrichment
            && signal_to_noise >= thresholds.min_signal_to_noise
            && duplicate_rate <= thresholds.max_duplicate_rate
            && final_reads >= thresholds.min_final_reads;
        
        if qc_pass {
            passed_qc_count += 1;
        }
        
        total_mapping_rate += mapping_rate;
        total_tss_enrichment += tss_enrichment;
        total_signal_to_noise += signal_to_noise;
        total_duplicate_rate += duplicate_rate;
        
        sample_metrics.push(SampleMetrics {
            sample_name: sample_name.clone(),
            mapping_rate,
            tss_enrichment,
            signal_to_noise,
            duplicate_rate,
            final_reads,
            quality_grade,
            qc_pass,
        });
    }
    
    let sample_count = sample_metrics.len();
    let successful_samples = sample_count; // All samples that have at least some data
    
    let overall_metrics = OverallMetrics {
        avg_mapping_rate: if sample_count > 0 { total_mapping_rate / sample_count as f64 } else { 0.0 },
        avg_tss_enrichment: if sample_count > 0 { total_tss_enrichment / sample_count as f64 } else { 0.0 },
        avg_signal_to_noise: if sample_count > 0 { total_signal_to_noise / sample_count as f64 } else { 0.0 },
        avg_duplicate_rate: if sample_count > 0 { total_duplicate_rate / sample_count as f64 } else { 0.0 },
        samples_passed_qc: passed_qc_count,
    };
    
    let recommendations = generate_recommendations(&overall_metrics, &sample_metrics, &thresholds);
    
    Ok(QCSummary {
        pipeline_info: PipelineInfo {
            version: pipeline_version.to_string(),
            analysis_date: analysis_date.to_string(),
            total_samples: sample_count,
            successful_samples,
        },
        sample_count,
        overall_metrics,
        sample_metrics,
        quality_thresholds: thresholds,
        recommendations,
    })
}

fn generate_recommendations(
    overall: &OverallMetrics,
    samples: &[SampleMetrics],
    thresholds: &QualityThresholds,
) -> Vec<String> {
    let mut recommendations = Vec::new();
    
    if overall.avg_mapping_rate < thresholds.min_mapping_rate {
        recommendations.push(format!(
            "Low average mapping rate ({:.1}%). Consider checking reference genome quality or sequencing adapter contamination.",
            overall.avg_mapping_rate
        ));
    }
    
    if overall.avg_tss_enrichment < thresholds.min_tss_enrichment {
        recommendations.push(format!(
            "Low TSS enrichment ({:.2}). This may indicate poor chromatin accessibility or library preparation issues.",
            overall.avg_tss_enrichment
        ));
    }
    
    if overall.avg_duplicate_rate > thresholds.max_duplicate_rate {
        recommendations.push(format!(
            "High duplicate rate ({:.1}%). Consider optimizing PCR amplification or increasing library complexity.",
            overall.avg_duplicate_rate
        ));
    }
    
    let failed_samples: Vec<_> = samples.iter().filter(|s| !s.qc_pass).collect();
    if !failed_samples.is_empty() {
        recommendations.push(format!(
            "{} samples failed QC criteria. Review individual sample metrics for detailed analysis.",
            failed_samples.len()
        ));
    }
    
    if overall.samples_passed_qc == samples.len() {
        recommendations.push("All samples passed QC criteria. Excellent data quality!".to_string());
    }
    
    if overall.avg_signal_to_noise < thresholds.min_signal_to_noise {
        recommendations.push(format!(
            "Low signal-to-noise ratio ({:.2}). Consider reviewing sequencing depth and library preparation protocols.",
            overall.avg_signal_to_noise
        ));
    }
    
    recommendations
}

fn write_metrics_table(
    summary: &QCSummary,
    output_file: &str,
    log_writer: &mut Option<BufWriter<File>>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut writer = BufWriter::new(File::create(output_file)?);
    
    // Write header
    writeln!(
        writer,
        "sample\tmapping_rate\ttss_enrichment\tsignal_to_noise\tduplicate_rate\tfinal_reads\tquality_grade\tqc_pass"
    )?;
    
    // Write sample data
    for sample in &summary.sample_metrics {
        writeln!(
            writer,
            "{}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{}\t{}\t{}",
            sample.sample_name,
            sample.mapping_rate,
            sample.tss_enrichment,
            sample.signal_to_noise,
            sample.duplicate_rate,
            sample.final_reads,
            sample.quality_grade,
            if sample.qc_pass { "PASS" } else { "FAIL" }
        )?;
    }
    
    writer.flush()?;
    log_message(log_writer, &format!("Metrics table written to: {}", output_file))?;
    
    Ok(())
}

fn write_html_report(
    summary: &QCSummary,
    output_file: &str,
    log_writer: &mut Option<BufWriter<File>>,
) -> Result<(), Box<dyn std::error::Error>> {
    let html_content = generate_html_report(summary)?;
    std::fs::write(output_file, html_content)?;
    log_message(log_writer, &format!("HTML report written to: {}", output_file))?;
    Ok(())
}

fn generate_html_report(summary: &QCSummary) -> Result<String, Box<dyn std::error::Error>> {
    let html = format!(r#"
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>ATAC-Seq QC Report</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
            color: #333;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            padding: 30px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        .header {{
            border-bottom: 3px solid #4CAF50;
            padding-bottom: 20px;
            margin-bottom: 30px;
        }}
        h1 {{
            color: #2E7D32;
            margin: 0;
            font-size: 2.5em;
        }}
        .subtitle {{
            color: #666;
            margin-top: 10px;
            font-size: 1.1em;
        }}
        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        .metric-card {{
            background: #f8f9fa;
            padding: 20px;
            border-radius: 8px;
            border-left: 4px solid #4CAF50;
        }}
        .metric-value {{
            font-size: 2em;
            font-weight: bold;
            color: #2E7D32;
        }}
        .metric-label {{
            color: #666;
            margin-top: 5px;
        }}
        .table-container {{
            overflow-x: auto;
            margin: 20px 0;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }}
        th, td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}
        th {{
            background-color: #4CAF50;
            color: white;
            font-weight: bold;
        }}
        tr:nth-child(even) {{
            background-color: #f2f2f2;
        }}
        tr:hover {{
            background-color: #e8f5e8;
        }}
        .pass {{
            color: #4CAF50;
            font-weight: bold;
        }}
        .fail {{
            color: #f44336;
            font-weight: bold;
        }}
        .recommendations {{
            background: #e3f2fd;
            padding: 20px;
            border-radius: 8px;
            border-left: 4px solid #2196F3;
            margin: 20px 0;
        }}
        .recommendations h3 {{
            color: #1976D2;
            margin-top: 0;
        }}
        .recommendations ul {{
            margin: 10px 0;
            padding-left: 20px;
        }}
        .recommendations li {{
            margin: 8px 0;
            line-height: 1.5;
        }}
        .quality-excellent {{ background-color: #e8f5e8; }}
        .quality-good {{ background-color: #fff3e0; }}
        .quality-acceptable {{ background-color: #fff8e1; }}
        .quality-poor {{ background-color: #ffebee; }}
        .quality-failed {{ background-color: #ffcdd2; }}
        .footer {{
            margin-top: 40px;
            padding-top: 20px;
            border-top: 1px solid #ddd;
            color: #666;
            text-align: center;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>ATAC-Seq Quality Control Report</h1>
            <div class="subtitle">
                Pipeline Version: {pipeline_version} | Analysis Date: {analysis_date}
            </div>
        </div>

        <div class="summary-grid">
            <div class="metric-card">
                <div class="metric-value">{total_samples}</div>
                <div class="metric-label">Total Samples</div>
            </div>
            <div class="metric-card">
                <div class="metric-value">{passed_samples}</div>
                <div class="metric-label">Samples Passed QC</div>
            </div>
            <div class="metric-card">
                <div class="metric-value">{avg_mapping_rate:.1}%</div>
                <div class="metric-label">Avg Mapping Rate</div>
            </div>
            <div class="metric-card">
                <div class="metric-value">{avg_tss_enrichment:.2}</div>
                <div class="metric-label">Avg TSS Enrichment</div>
            </div>
        </div>

        <h2>Quality Thresholds</h2>
        <div class="table-container">
            <table>
                <thead>
                    <tr>
                        <th>Metric</th>
                        <th>Threshold</th>
                        <th>Current Average</th>
                        <th>Status</th>
                    </tr>
                </thead>
                <tbody>
                    <tr>
                        <td>Mapping Rate</td>
                        <td>≥ {min_mapping_rate}%</td>
                        <td>{avg_mapping_rate:.1}%</td>
                        <td class="{mapping_status}">{mapping_status_text}</td>
                    </tr>
                    <tr>
                        <td>TSS Enrichment</td>
                        <td>≥ {min_tss_enrichment}</td>
                        <td>{avg_tss_enrichment:.2}</td>
                        <td class="{tss_status}">{tss_status_text}</td>
                    </tr>
                    <tr>
                        <td>Signal-to-Noise</td>
                        <td>≥ {min_signal_to_noise}</td>
                        <td>{avg_signal_to_noise:.2}</td>
                        <td class="{snr_status}">{snr_status_text}</td>
                    </tr>
                    <tr>
                        <td>Duplicate Rate</td>
                        <td>≤ {max_duplicate_rate}%</td>
                        <td>{avg_duplicate_rate:.1}%</td>
                        <td class="{dup_status}">{dup_status_text}</td>
                    </tr>
                </tbody>
            </table>
        </div>

        <h2>Sample Metrics</h2>
        <div class="table-container">
            <table>
                <thead>
                    <tr>
                        <th>Sample</th>
                        <th>Mapping Rate (%)</th>
                        <th>TSS Enrichment</th>
                        <th>Signal-to-Noise</th>
                        <th>Duplicate Rate (%)</th>
                        <th>Final Reads</th>
                        <th>Quality Grade</th>
                        <th>QC Status</th>
                    </tr>
                </thead>
                <tbody>
                    {sample_rows}
                </tbody>
            </table>
        </div>

        <div class="recommendations">
            <h3>Recommendations</h3>
            <ul>
                {recommendations_list}
            </ul>
        </div>

        <div class="footer">
            <p>Generated by ATAC-Seq QC Pipeline v{pipeline_version}</p>
            <p>For questions or support, please contact your bioinformatics team.</p>
        </div>
    </div>
</body>
</html>
"#,
        pipeline_version = summary.pipeline_info.version,
        analysis_date = summary.pipeline_info.analysis_date,
        total_samples = summary.sample_count,
        passed_samples = summary.overall_metrics.samples_passed_qc,
        avg_mapping_rate = summary.overall_metrics.avg_mapping_rate,
        avg_tss_enrichment = summary.overall_metrics.avg_tss_enrichment,
        avg_signal_to_noise = summary.overall_metrics.avg_signal_to_noise,
        avg_duplicate_rate = summary.overall_metrics.avg_duplicate_rate,
        min_mapping_rate = summary.quality_thresholds.min_mapping_rate,
        min_tss_enrichment = summary.quality_thresholds.min_tss_enrichment,
        min_signal_to_noise = summary.quality_thresholds.min_signal_to_noise,
        max_duplicate_rate = summary.quality_thresholds.max_duplicate_rate,
        mapping_status = if summary.overall_metrics.avg_mapping_rate >= summary.quality_thresholds.min_mapping_rate { "pass" } else { "fail" },
        mapping_status_text = if summary.overall_metrics.avg_mapping_rate >= summary.quality_thresholds.min_mapping_rate { "PASS" } else { "FAIL" },
        tss_status = if summary.overall_metrics.avg_tss_enrichment >= summary.quality_thresholds.min_tss_enrichment { "pass" } else { "fail" },
        tss_status_text = if summary.overall_metrics.avg_tss_enrichment >= summary.quality_thresholds.min_tss_enrichment { "PASS" } else { "FAIL" },
        snr_status = if summary.overall_metrics.avg_signal_to_noise >= summary.quality_thresholds.min_signal_to_noise { "pass" } else { "fail" },
        snr_status_text = if summary.overall_metrics.avg_signal_to_noise >= summary.quality_thresholds.min_signal_to_noise { "PASS" } else { "FAIL" },
        dup_status = if summary.overall_metrics.avg_duplicate_rate <= summary.quality_thresholds.max_duplicate_rate { "pass" } else { "fail" },
        dup_status_text = if summary.overall_metrics.avg_duplicate_rate <= summary.quality_thresholds.max_duplicate_rate { "PASS" } else { "FAIL" },
        sample_rows = generate_sample_rows(&summary.sample_metrics),
        recommendations_list = generate_recommendations_html(&summary.recommendations),
    );

    Ok(html)
}

fn generate_sample_rows(samples: &[SampleMetrics]) -> String {
    samples
        .iter()
        .map(|sample| {
            let quality_class = match sample.quality_grade.as_str() {
                "Excellent" => "quality-excellent",
                "Good" => "quality-good", 
                "Acceptable" => "quality-acceptable",
                "Poor" => "quality-poor",
                _ => "quality-failed",
            };
            
            let qc_class = if sample.qc_pass { "pass" } else { "fail" };
            let qc_text = if sample.qc_pass { "PASS" } else { "FAIL" };
            
            format!(
                r#"<tr class="{quality_class}">
                    <td><strong>{sample_name}</strong></td>
                    <td>{mapping_rate:.2}</td>
                    <td>{tss_enrichment:.2}</td>
                    <td>{signal_to_noise:.2}</td>
                    <td>{duplicate_rate:.2}</td>
                    <td>{final_reads:,}</td>
                    <td>{quality_grade}</td>
                    <td class="{qc_class}"><strong>{qc_text}</strong></td>
                </tr>"#,
                quality_class = quality_class,
                sample_name = sample.sample_name,
                mapping_rate = sample.mapping_rate,
                tss_enrichment = sample.tss_enrichment,
                signal_to_noise = sample.signal_to_noise,
                duplicate_rate = sample.duplicate_rate,
                final_reads = sample.final_reads,
                quality_grade = sample.quality_grade,
                qc_class = qc_class,
                qc_text = qc_text,
            )
        })
        .collect::<Vec<_>>()
        .join("\n")
}

fn generate_recommendations_html(recommendations: &[String]) -> String {
    recommendations
        .iter()
        .map(|rec| format!("<li>{}</li>", rec))
        .collect::<Vec<_>>()
        .join("\n")
}