#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Fragment Shift Estimation and Peak Calling Pipeline
params.input_dir = './data'
params.output_dir = './results'
params.reference_genome = './reference/genome.fa'

// Shift estimation parameters
params.max_shift = 500
params.min_mapq = 10
params.bin_size = 1
params.use_fft = true
params.smoothing_window = 5
params.sampling_factor = 1

// Peak calling parameters
params.window_size = 200
params.min_coverage = 5.0
params.pvalue_threshold = 0.05
params.min_peak_width = 50
params.max_peak_width = 5000

// Performance parameters
params.threads = 4
params.memory = '8 GB'
params.time_limit = '4h'

// Analysis scope (optional chromosome-specific analysis)
params.target_chromosome = null
params.analysis_start = null
params.analysis_end = null

log.info """
===============================================
Fragment Shift Estimation and Peak Calling Pipeline
===============================================
Input directory     : ${params.input_dir}
Output directory    : ${params.output_dir}
Reference genome    : ${params.reference_genome}

Shift Estimation:
- Max shift range   : ${params.max_shift} bp
- Bin size         : ${params.bin_size} bp
- Use FFT          : ${params.use_fft}
- Min MAPQ         : ${params.min_mapq}
- Target chromosome: ${params.target_chromosome ?: 'All chromosomes'}

Peak Calling:
- Window size      : ${params.window_size} bp
- Min coverage     : ${params.min_coverage}
- P-value threshold: ${params.pvalue_threshold}
- Peak width range : ${params.min_peak_width}-${params.max_peak_width} bp

Performance:
- Threads          : ${params.threads}
- Memory           : ${params.memory}
===============================================
"""

// Validate inputs
process VALIDATE_SHIFT_INPUTS {
    tag "shift_validation"
    
    input:
    path input_dir
    path reference_genome
    
    output:
    stdout
    
    script:
    """
    echo "=== Shift Pipeline Input Validation ==="
    
    # Check input directory
    test -d "${input_dir}" || (echo "ERROR: Input directory not found: ${input_dir}" && exit 1)
    
    # Check reference genome
    test -f "${reference_genome}" || (echo "ERROR: Reference genome not found: ${reference_genome}" && exit 1)
    
    # Count FASTQ files
    fastq_count=\$(find ${input_dir} -name "*.fastq.gz" -o -name "*.fq.gz" | wc -l)
    test \$fastq_count -gt 0 || (echo "ERROR: No FASTQ files found in ${input_dir}" && exit 1)
    
    echo "✓ Found \$fastq_count FASTQ files"
    echo "✓ All inputs validated successfully"
    """
}

// Build reference indices
process BUILD_SHIFT_INDICES {
    tag "shift_indexing"
    publishDir "${params.output_dir}/reference", mode: 'copy'
    
    input:
    path reference_genome
    
    output:
    path "genome_index*", emit: bt2_index
    path "${reference_genome}.fai", emit: fai_index
    
    script:
    """
    echo "Building reference indices for shift analysis..."
    
    # Build bowtie2 index
    bowtie2-build --threads ${params.threads} ${reference_genome} genome_index
    
    # Build samtools index
    samtools faidx ${reference_genome}
    
    echo "Reference indexing completed"
    """
}

// Quality control and trimming
process FASTQ_PREPROCESSING {
    tag "${sample_id}"
    publishDir "${params.output_dir}/qc", mode: 'copy', pattern: "*.{json,html}"
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq.gz"), emit: trimmed_reads
    tuple val(sample_id), path("${sample_id}_fastp.{json,html}"), emit: qc_reports
    
    script:
    """
    # Enhanced preprocessing for shift analysis
    fastp \\
        -i ${reads} \\
        -o ${sample_id}_trimmed.fastq.gz \\
        --qualified_quality_phred 20 \\
        --unqualified_percent_limit 20 \\
        --length_required 36 \\
        --detect_adapter_for_pe \\
        --correction \\
        --thread ${params.threads} \\
        --json ${sample_id}_fastp.json \\
        --html ${sample_id}_fastp.html
    """
}

// Optimized alignment for cross-correlation analysis
process SHIFT_ALIGNMENT {
    tag "${sample_id}"
    publishDir "${params.output_dir}/alignments", mode: 'copy', pattern: "*.{bam,bai}"
    
    input:
    tuple val(sample_id), path(trimmed_reads)
    path bt2_index
    
    output:
    tuple val(sample_id), path("${sample_id}_aligned.bam"), emit: aligned_bam
    tuple val(sample_id), path("${sample_id}_aligned.bam.bai"), emit: bam_index
    tuple val(sample_id), path("${sample_id}_alignment_stats.txt"), emit: alignment_stats
    
    script:
    """
    echo "Aligning ${sample_id} for cross-correlation analysis..."
    
    # Alignment optimized for fragment shift estimation
    bowtie2 \\
        -x genome_index \\
        -U ${trimmed_reads} \\
        --threads ${params.threads} \\
        --very-sensitive \\
        --no-discordant \\
        --no-mixed \\
        -k 1 \\
        2> ${sample_id}_bowtie2.log \\
        | samtools view -bS -F 4 -q ${params.min_mapq} - \\
        | samtools sort -@ ${params.threads} -o ${sample_id}_aligned.bam -
    
    # Index BAM file
    samtools index ${sample_id}_aligned.bam
    
    # Generate alignment statistics
    samtools flagstat ${sample_id}_aligned.bam > ${sample_id}_alignment_stats.txt
    samtools stats ${sample_id}_aligned.bam >> ${sample_id}_alignment_stats.txt
    
    # Extract key metrics for shift analysis
    echo "=== Alignment Quality for Shift Analysis ===" >> ${sample_id}_alignment_stats.txt
    echo "Total aligned reads: \$(samtools view -c -F 4 ${sample_id}_aligned.bam)" >> ${sample_id}_alignment_stats.txt
    echo "Forward strand reads: \$(samtools view -c -F 20 ${sample_id}_aligned.bam)" >> ${sample_id}_alignment_stats.txt
    echo "Reverse strand reads: \$(samtools view -c -f 16 -F 4 ${sample_id}_aligned.bam)" >> ${sample_id}_alignment_stats.txt
    echo "High quality reads (MAPQ>=${params.min_mapq}): \$(samtools view -c -F 4 -q ${params.min_mapq} ${sample_id}_aligned.bam)" >> ${sample_id}_alignment_stats.txt
    """
}

// Fragment shift estimation using cross-correlation
process ESTIMATE_FRAGMENT_SHIFT {
    tag "${sample_id}"
    publishDir "${params.output_dir}/shift_analysis", mode: 'copy'
    
    input:
    tuple val(sample_id), path(aligned_bam)
    
    output:
    tuple val(sample_id), path("${sample_id}_shift_estimate.json"), emit: shift_json
    tuple val(sample_id), path("${sample_id}_correlation_profile.json"), emit: correlation_profile
    tuple val(sample_id), path("${sample_id}_shift_analysis.log"), emit: shift_log
    
    script:
    def chrom_option = params.target_chromosome ? "--chromosome ${params.target_chromosome}" : ""
    def start_option = params.analysis_start ? "--start-pos ${params.analysis_start}" : ""
    def end_option = params.analysis_end ? "--end-pos ${params.analysis_end}" : ""
    def fft_option = params.use_fft ? "--use-fft" : ""
    
    """
    echo "Estimating fragment shift for ${sample_id}..."
    
    # Run cross-correlation analysis
    cargo run --release --bin shift_estimator -- \\
        --input ${aligned_bam} \\
        --output ${sample_id}_shift_estimate.json \\
        --correlation-output ${sample_id}_correlation_profile.json \\
        --max-shift ${params.max_shift} \\
        --min-mapq ${params.min_mapq} \\
        --bin-size ${params.bin_size} \\
        --smoothing-window ${params.smoothing_window} \\
        --sampling-factor ${params.sampling_factor} \\
        ${chrom_option} \\
        ${start_option} \\
        ${end_option} \\
        ${fft_option} \\
        --threads ${params.threads} \\
        --verbose \\
        2>&1 | tee ${sample_id}_shift_analysis.log
    
    # Validate shift estimation
    if [ ! -f "${sample_id}_shift_estimate.json" ]; then
        echo "ERROR: Shift estimation failed for ${sample_id}"
        exit 1
    fi
    
    # Extract key metrics for reporting
    echo "" >> ${sample_id}_shift_analysis.log
    echo "=== Shift Estimation Summary ===" >> ${sample_id}_shift_analysis.log
    echo "Sample: ${sample_id}" >> ${sample_id}_shift_analysis.log
    echo "Estimated shift: \$(jq -r '.estimated_shift' ${sample_id}_shift_estimate.json) bp" >> ${sample_id}_shift_analysis.log
    echo "Confidence score: \$(jq -r '.confidence_score' ${sample_id}_shift_estimate.json)" >> ${sample_id}_shift_analysis.log
    echo "Signal-to-noise ratio: \$(jq -r '.signal_to_noise_ratio' ${sample_id}_shift_estimate.json)" >> ${sample_id}_shift_analysis.log
    """
}

// Shift-aware peak calling
process SHIFT_AWARE_PEAK_CALLING {
    tag "${sample_id}"
    publishDir "${params.output_dir}/peaks", mode: 'copy'
    
    input:
    tuple val(sample_id), path(aligned_bam), path(shift_json)
    
    output:
    tuple val(sample_id), path("${sample_id}_peaks.bed"), emit: peaks_bed
    tuple val(sample_id), path("${sample_id}_peak_calling.log"), emit: peak_log
    tuple val(sample_id), path("${sample_id}_peak_summary.txt"), emit: peak_summary
    
    script:
    """
    echo "Calling peaks with shift correction for ${sample_id}..."
    
    # Run shift-aware peak calling
    cargo run --release --bin peak_caller -- \\
        ${aligned_bam} \\
        ${shift_json} \\
        --output ${sample_id}_peaks.bed \\
        --window-size ${params.window_size} \\
        --min-coverage ${params.min_coverage} \\
        --pvalue-threshold ${params.pvalue_threshold} \\
        --min-peak-width ${params.min_peak_width} \\
        --max-peak-width ${params.max_peak_width} \\
        --min-mapq ${params.min_mapq} \\
        --threads ${params.threads} \\
        --verbose \\
        2>&1 | tee ${sample_id}_peak_calling.log
    
    # Generate peak summary
    echo "Peak Calling Summary for ${sample_id}" > ${sample_id}_peak_summary.txt
    echo "=====================================" >> ${sample_id}_peak_summary.txt
    echo "Date: \$(date)" >> ${sample_id}_peak_summary.txt
    echo "" >> ${sample_id}_peak_summary.txt
    
    # Count peaks
    peak_count=\$(grep -v "^#" ${sample_id}_peaks.bed | wc -l)
    echo "Total peaks called: \$peak_count" >> ${sample_id}_peak_summary.txt
    
    if [ \$peak_count -gt 0 ]; then
        # Calculate peak statistics
        echo "Peak width statistics:" >> ${sample_id}_peak_summary.txt
        awk 'NR>2 {print \$3-\$2}' ${sample_id}_peaks.bed | awk '
            {
                widths[NR] = \$1
                sum += \$1
                if (NR == 1) {min = max = \$1}
                if (\$1 < min) min = \$1
                if (\$1 > max) max = \$1
            }
            END {
                if (NR > 0) {
                    mean = sum / NR
                    print "  Mean width: " mean " bp"
                    print "  Min width: " min " bp"
                    print "  Max width: " max " bp"
                }
            }
        ' >> ${sample_id}_peak_summary.txt
        
        # Signal statistics
        echo "" >> ${sample_id}_peak_summary.txt
        echo "Signal statistics:" >> ${sample_id}_peak_summary.txt
        awk 'NR>2 {print \$7}' ${sample_id}_peaks.bed | awk '
            {
                sum += \$1
                if (NR == 1) {min = max = \$1}
                if (\$1 < min) min = \$1
                if (\$1 > max) max = \$1
            }
            END {
                if (NR > 0) {
                    mean = sum / NR
                    print "  Mean signal: " mean
                    print "  Min signal: " min
                    print "  Max signal: " max
                }
            }
        ' >> ${sample_id}_peak_summary.txt
        
        # Significance statistics
        echo "" >> ${sample_id}_peak_summary.txt
        echo "Significance statistics:" >> ${sample_id}_peak_summary.txt
        echo "  Peaks with p-value < 0.001: \$(awk 'NR>2 && \$8 < 0.001' ${sample_id}_peaks.bed | wc -l)" >> ${sample_id}_peak_summary.txt
        echo "  Peaks with p-value < 0.01: \$(awk 'NR>2 && \$8 < 0.01' ${sample_id}_peaks.bed | wc -l)" >> ${sample_id}_peak_summary.txt
        echo "  Peaks with p-value < 0.05: \$(awk 'NR>2 && \$8 < 0.05' ${sample_id}_peaks.bed | wc -l)" >> ${sample_id}_peak_summary.txt
    else
        echo "No peaks found - check data quality and parameters" >> ${sample_id}_peak_summary.txt
    fi
    
    # Add shift information
    echo "" >> ${sample_id}_peak_summary.txt
    echo "Fragment shift information:" >> ${sample_id}_peak_summary.txt
    echo "  Applied shift: \$(jq -r '.estimated_shift' ${shift_json}) bp" >> ${sample_id}_peak_summary.txt
    echo "  Shift confidence: \$(jq -r '.confidence_score' ${shift_json})" >> ${sample_id}_peak_summary.txt
    """
}

// Comparative analysis across samples
process COMPARE_SHIFT_ESTIMATES {
    tag "shift_comparison"
    publishDir "${params.output_dir}/analysis", mode: 'copy'
    
    input:
    path shift_jsons
    path correlation_profiles
    
    output:
    path "shift_comparison.txt", emit: shift_comparison
    path "shift_analysis_summary.json", emit: shift_summary
    
    script:
    """
    echo "Comparing fragment shift estimates across samples..."
    
    python3 << 'EOF'
import json
import os
import statistics
from datetime import datetime

# Collect shift estimates
shift_data = []
for json_file in [f for f in os.listdir('.') if f.endswith('_shift_estimate.json')]:
    sample_id = json_file.replace('_shift_estimate.json', '')
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
            data['sample_id'] = sample_id
            shift_data.append(data)
    except Exception as e:
        print(f"Error loading {json_file}: {e}")

if not shift_data:
    print("No shift data found")
    with open('shift_comparison.txt', 'w') as f:
        f.write('No shift estimates available\\n')
    with open('shift_analysis_summary.json', 'w') as f:
        json.dump({'error': 'No data'}, f)
    exit(0)

# Generate comparison report
shifts = [d['estimated_shift'] for d in shift_data]
confidences = [d['confidence_score'] for d in shift_data]
snrs = [d['signal_to_noise_ratio'] for d in shift_data]

comparison_text = f"""Fragment Shift Analysis Summary
Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
Total samples analyzed: {len(shift_data)}

Fragment Shift Statistics:
- Mean shift: {statistics.mean(shifts):.1f} bp
- Median shift: {statistics.median(shifts):.1f} bp
- Standard deviation: {statistics.stdev(shifts) if len(shifts) > 1 else 0:.1f} bp
- Range: {min(shifts):.0f} - {max(shifts):.0f} bp

Confidence Score Statistics:
- Mean confidence: {statistics.mean(confidences):.3f}
- Median confidence: {statistics.median(confidences):.3f}
- Samples with high confidence (>0.7): {sum(1 for c in confidences if c > 0.7)}
- Samples with low confidence (<0.5): {sum(1 for c in confidences if c < 0.5)}

Signal-to-Noise Ratio Statistics:
- Mean SNR: {statistics.mean(snrs):.2f}
- Median SNR: {statistics.median(snrs):.2f}
- Samples with good SNR (>3.0): {sum(1 for s in snrs if s > 3.0)}

Individual Sample Results:
{'='*50}
"""

for data in sorted(shift_data, key=lambda x: x['sample_id']):
    comparison_text += f"""
Sample: {data['sample_id']}
  Estimated shift: {data['estimated_shift']} bp
  Confidence score: {data['confidence_score']:.3f}
  Signal-to-noise ratio: {data['signal_to_noise_ratio']:.2f}
  Correlation peak: {data['correlation_peak']:.4f}
  Processing time: {data['processing_stats']['analysis_time_seconds']:.1f}s
"""

# Quality assessment
quality_issues = []
for data in shift_data:
    sample = data['sample_id']
    if data['confidence_score'] < 0.5:
        quality_issues.append(f"  - {sample}: Low confidence ({data['confidence_score']:.3f})")
    if data['signal_to_noise_ratio'] < 2.0:
        quality_issues.append(f"  - {sample}: Low SNR ({data['signal_to_noise_ratio']:.2f})")

if quality_issues:
    comparison_text += f"""

Quality Concerns:
{chr(10).join(quality_issues)}

Recommendations:
- Check data quality for samples with low confidence/SNR
- Consider increasing sequencing depth
- Verify ChIP efficiency and antibody specificity
"""
else:
    comparison_text += """

✓ All samples show good shift estimation quality
"""

# Write comparison report
with open('shift_comparison.txt', 'w') as f:
    f.write(comparison_text)

# Create summary JSON
summary = {
    'analysis_date': datetime.now().isoformat(),
    'total_samples': len(shift_data),
    'shift_statistics': {
        'mean': statistics.mean(shifts),
        'median': statistics.median(shifts),
        'std': statistics.stdev(shifts) if len(shifts) > 1 else 0,
        'min': min(shifts),
        'max': max(shifts)
    },
    'confidence_statistics': {
        'mean': statistics.mean(confidences),
        'median': statistics.median(confidences),
        'high_confidence_count': sum(1 for c in confidences if c > 0.7),
        'low_confidence_count': sum(1 for c in confidences if c < 0.5)
    },
    'snr_statistics': {
        'mean': statistics.mean(snrs),
        'median': statistics.median(snrs),
        'good_snr_count': sum(1 for s in snrs if s > 3.0)
    },
    'sample_details': shift_data
}

with open('shift_analysis_summary.json', 'w') as f:
    json.dump(summary, f, indent=2)

print("Shift comparison analysis completed")
EOF
    """
}

// Aggregate peak calling results
process AGGREGATE_PEAK_RESULTS {
    tag "peak_aggregation"
    publishDir "${params.output_dir}/summary", mode: 'copy'
    
    input:
    path peak_beds
    path peak_summaries
    path shift_comparison
    
    output:
    path "peak_calling_report.html", emit: peak_report
    path "merged_peaks.bed", emit: merged_peaks
    path "pipeline_summary.txt", emit: pipeline_summary
    
    script:
    """
    echo "Aggregating peak calling results..."
    
    python3 << 'EOF'
import os
import json
from datetime import datetime

# Collect peak statistics
peak_stats = []
total_peaks = 0

for summary_file in [f for f in os.listdir('.') if f.endswith('_peak_summary.txt')]:
    sample_id = summary_file.replace('_peak_summary.txt', '')
    try:
        with open(summary_file, 'r') as f:
            content = f.read()
            # Extract peak count
            for line in content.split('\\n'):
                if 'Total peaks called:' in line:
                    count = int(line.split(':')[1].strip())
                    peak_stats.append((sample_id, count))
                    total_peaks += count
                    break
    except Exception as e:
        print(f"Error processing {summary_file}: {e}")

# Generate HTML report
html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>Fragment Shift and Peak Calling Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; background: #f8f9fa; }}
        .container {{ max-width: 1200px; margin: 0 auto; background: white; padding: 30px; border-radius: 10px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); }}
        h1, h2 {{ color: #2c3e50; }}
        h1 {{ border-bottom: 3px solid #3498db; padding-bottom: 10px; }}
        .summary-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 20px; margin: 20px 0; }}
        .metric-card {{ background: #ecf0f1; padding: 20px; border-radius: 8px; border-left: 5px solid #3498db; }}
        .metric-value {{ font-size: 1.8em; font-weight: bold; color: #2c3e50; }}
        .metric-label {{ color: #7f8c8d; font-size: 0.9em; }}
        table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 12px; text-align: left; }}
        th {{ background: #34495e; color: white; }}
        tr:nth-child(even) {{ background: #f2f2f2; }}
        .high-quality {{ color: #27ae60; font-weight: bold; }}
        .medium-quality {{ color: #f39c12; font-weight: bold; }}
        .low-quality {{ color: #e74c3c; font-weight: bold; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 Fragment Shift and Peak Calling Analysis Report</h1>
        
        <div class="summary-grid">
            <div class="metric-card">
                <div class="metric-value">{len(peak_stats)}</div>
                <div class="metric-label">Samples Analyzed</div>
            </div>
            <div class="metric-card">
                <div class="metric-value">{total_peaks:,}</div>
                <div class="metric-label">Total Peaks Called</div>
            </div>
            <div class="metric-card">
                <div class="metric-value">{total_peaks // len(peak_stats) if peak_stats else 0:,}</div>
                <div class="metric-label">Average Peaks per Sample</div>
            </div>
            <div class="metric-card">
                <div class="metric-value">{datetime.now().strftime('%Y-%m-%d')}</div>
                <div class="metric-label">Analysis Date</div>
            </div>
        </div>

        <h2>📊 Sample Results</h2>
        <table>
            <tr>
                <th>Sample ID</th>
                <th>Peaks Called</th>
                <th>Quality Assessment</th>
            </tr>
"""

# Add sample rows
for sample_id, peak_count in sorted(peak_stats):
    quality_class = "high-quality" if peak_count > 1000 else "medium-quality" if peak_count > 100 else "low-quality"
    quality_text = "High" if peak_count > 1000 else "Medium" if peak_count > 100 else "Low"
    
    html_content += f"""
            <tr>
                <td>{sample_id}</td>
                <td>{peak_count:,}</td>
                <td><span class="{quality_class}">{quality_text}</span></td>
            </tr>
    """

# Load shift comparison if available
shift_summary = ""
try:
    with open('shift_comparison.txt', 'r') as f:
        shift_summary = f.read().replace('\\n', '<br>')
except:
    shift_summary = "Shift comparison data not available"

html_content += f"""
        </table>

        <h2>🔬 Fragment Shift Analysis</h2>
        <div class="metric-card">
            <pre style="white-space: pre-wrap; font-family: monospace; font-size: 0.9em;">
{shift_summary.replace('<br>', chr(10))}
            </pre>
        </div>

        <h2>📋 Analysis Parameters</h2>
        <div class="metric-card">
            <strong>Shift Estimation:</strong><br>
            • Max shift range: {params.max_shift} bp<br>
            • Bin size: {params.bin_size} bp<br>
            • FFT acceleration: {params.use_fft}<br>
            • Min MAPQ: {params.min_mapq}<br><br>
            
            <strong>Peak Calling:</strong><br>
            • Window size: {params.window_size} bp<br>
            • Min coverage: {params.min_coverage}<br>
            • P-value threshold: {params.pvalue_threshold}<br>
            • Peak width: {params.min_peak_width}-{params.max_peak_width} bp
        </div>

        <footer style="margin-top: 40px; padding-top: 20px; border-top: 2px solid #ecf0f1; text-align: center; color: #7f8c8d;">
            <p>Generated by Fragment Shift Estimation Pipeline | {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </footer>
    </div>
</body>
</html>
"""

# Write HTML report
with open('peak_calling_report.html', 'w') as f:
    f.write(html_content)

print("Peak calling report generated")
EOF

    # Merge all peak BED files
    echo "Merging peak files..."
    echo "# Merged peaks from all samples" > merged_peaks.bed
    echo "# Generated: \$(date)" >> merged_peaks.bed
    echo "# chrom\tstart\tend\tname\tscore\tstrand\tsignal\tpValue\tqValue\tsummit" >> merged_peaks.bed
    
    for bed_file in *_peaks.bed; do
        if [ -f "\$bed_file" ]; then
            sample=\$(basename \$bed_file _peaks.bed)
            grep -v "^#" \$bed_file | awk -v sample="\$sample" 'BEGIN{OFS="\\t"} {print \$1,\$2,\$3,sample"_"\$4,\$5,\$6,\$7,\$8,\$9,\$10}' >> merged_peaks.bed
        fi
    done
    
    # Generate pipeline summary
    cat > pipeline_summary.txt << EOF
Fragment Shift Estimation and Peak Calling Pipeline Summary
==========================================================
Date: \$(date)
Total samples processed: ${peak_stats.length if 'peak_stats' in locals() else 0}
Total peaks called: ${total_peaks if 'total_peaks' in locals() else 0}

Pipeline Parameters:
- Max shift range: ${params.max_shift} bp
- Peak calling window: ${params.window_size} bp
- Significance threshold: ${params.pvalue_threshold}
- FFT acceleration: ${params.use_fft}

Output Files:
- Individual peaks: ${params.output_dir}/peaks/
- Shift estimates: ${params.output_dir}/shift_analysis/
- Merged peaks: ${params.output_dir}/summary/merged_peaks.bed
- Analysis report: ${params.output_dir}/summary/peak_calling_report.html

For detailed results, see the HTML report.
EOF
    """
}

// Main workflow
workflow {
    // Input validation
    VALIDATE_SHIFT_INPUTS(params.input_dir, params.reference_genome)
    
    // Create input channels
    fastq_ch = Channel
        .fromPath("${params.input_dir}/*.{fastq,fq}.gz")
        .map { file -> 
            def sample_id = file.baseName.toString().replaceAll(/\.(fastq|fq)$/, '')
            tuple(sample_id, file)
        }
    
    // Build reference indices
    BUILD_SHIFT_INDICES(params.reference_genome)
    
    // Preprocess FASTQ files
    FASTQ_PREPROCESSING(fastq_ch)
    
    // Alignment optimized for shift analysis
    SHIFT_ALIGNMENT(FASTQ_PREPROCESSING.out.trimmed_reads, BUILD_SHIFT_INDICES.out.bt2_index)
    
    // Fragment shift estimation
    ESTIMATE_FRAGMENT_SHIFT(SHIFT_ALIGNMENT.out.aligned_bam)
    
    // Shift-aware peak calling
    peak_input = SHIFT_ALIGNMENT.out.aligned_bam.join(ESTIMATE_FRAGMENT_SHIFT.out.shift_json)
    SHIFT_AWARE_PEAK_CALLING(peak_input)
    
    // Comparative analysis
    COMPARE_SHIFT_ESTIMATES(
        ESTIMATE_FRAGMENT_SHIFT.out.shift_json.collect(),
        ESTIMATE_FRAGMENT_SHIFT.out.correlation_profile.collect()
    )
    
    // Aggregate results
    AGGREGATE_PEAK_RESULTS(
        SHIFT_AWARE_PEAK_CALLING.out.peaks_bed.collect(),
        SHIFT_AWARE_PEAK_CALLING.out.peak_summary.collect(),
        COMPARE_SHIFT_ESTIMATES.out.shift_comparison
    )
}

workflow.onComplete {
    log.info """
    ================================================
    🎉 Fragment Shift Pipeline Complete!
    ================================================
    Success: ${workflow.success}
    Duration: ${workflow.duration}
    Start time: ${workflow.start}
    End time: ${workflow.complete}
    Work directory: ${workflow.workDir}
    Results directory: ${params.output_dir}
    
    📊 Key Outputs:
    - Shift estimates: ${params.output_dir}/shift_analysis/
    - Peak calls: ${params.output_dir}/peaks/
    - Analysis report: ${params.output_dir}/summary/peak_calling_report.html
    - Merged peaks: ${params.output_dir}/summary/merged_peaks.bed
    
    🔬 Analysis Summary:
    - Fragment shift estimation with cross-correlation
    - FFT acceleration: ${params.use_fft}
    - Shift-corrected peak calling
    - Multi-sample comparison and validation
    ================================================
    """.stripIndent()
}