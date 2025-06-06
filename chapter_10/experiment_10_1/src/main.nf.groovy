#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Pipeline parameters with defaults
params.input_dir = './data'
params.output_dir = './results'
params.reference_genome = './reference/genome.fa'
params.adapter_file = null
params.min_quality = 20
params.window_size = 200
params.min_coverage = 5.0
params.pvalue_threshold = 0.05
params.threads = 4

// Print pipeline parameters
log.info """
========================================
Epigenomic Peak Calling Pipeline
========================================
Input directory    : ${params.input_dir}
Output directory   : ${params.output_dir}
Reference genome   : ${params.reference_genome}
Adapter file       : ${params.adapter_file ?: 'Auto-detect'}
Min quality        : ${params.min_quality}
Window size        : ${params.window_size}
Min coverage       : ${params.min_coverage}
P-value threshold  : ${params.pvalue_threshold}
Threads            : ${params.threads}
========================================
"""

process VALIDATE_INPUTS {
    tag "validation"
    
    input:
    path input_dir
    path reference_genome
    
    output:
    stdout
    
    script:
    """
    echo "Validating inputs..."
    
    if [ ! -d "${input_dir}" ]; then
        echo "ERROR: Input directory does not exist: ${input_dir}"
        exit 1
    fi
    
    if [ ! -f "${reference_genome}" ]; then
        echo "ERROR: Reference genome file does not exist: ${reference_genome}"
        exit 1
    fi
    
    fastq_count=\$(find ${input_dir} -name "*.fastq.gz" -o -name "*.fq.gz" | wc -l)
    if [ \$fastq_count -eq 0 ]; then
        echo "ERROR: No FASTQ files found in ${input_dir}"
        exit 1
    fi
    
    echo "Validation passed. Found \$fastq_count FASTQ files."
    """
}

process QC_FASTQC {
    tag "${sample_id}"
    publishDir "${params.output_dir}/qc/fastqc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("*_fastqc.{html,zip}")
    
    script:
    """
    fastqc -t ${params.threads} -o . ${reads}
    """
}

process QC_AND_TRIM {
    tag "${sample_id}"
    publishDir "${params.output_dir}/qc/trimmed", mode: 'copy', pattern: "*.{json,html}"
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq.gz"), emit: trimmed_reads
    tuple val(sample_id), path("${sample_id}_fastp.{json,html}"), emit: qc_reports
    
    script:
    def adapter_option = params.adapter_file ? "--adapter_sequence_r1 \$(head -n 2 ${params.adapter_file} | tail -n 1)" : ""
    """
    fastp \\
        -i ${reads} \\
        -o ${sample_id}_trimmed.fastq.gz \\
        --qualified_quality_phred ${params.min_quality} \\
        --length_required 36 \\
        --thread ${params.threads} \\
        --json ${sample_id}_fastp.json \\
        --html ${sample_id}_fastp.html \\
        ${adapter_option}
    """
}

process BUILD_INDEX {
    tag "reference_indexing"
    publishDir "${params.output_dir}/reference", mode: 'copy'
    
    input:
    path reference_genome
    
    output:
    path "genome_index*", emit: bt2_index
    path "${reference_genome}.fai", emit: fai_index
    
    script:
    """
    # Build bowtie2 index
    bowtie2-build --threads ${params.threads} ${reference_genome} genome_index
    
    # Build samtools index
    samtools faidx ${reference_genome}
    """
}

process ALIGN {
    tag "${sample_id}"
    publishDir "${params.output_dir}/alignments", mode: 'copy', pattern: "*.{bam,bai,stats}"
    
    input:
    tuple val(sample_id), path(trimmed_reads)
    path bt2_index
    
    output:
    tuple val(sample_id), path("${sample_id}_aligned.bam"), emit: aligned_bam
    tuple val(sample_id), path("${sample_id}_aligned.bam.bai"), emit: bam_index
    tuple val(sample_id), path("${sample_id}_alignment_stats.txt"), emit: stats
    
    script:
    """
    # Align reads
    bowtie2 \\
        -x genome_index \\
        -U ${trimmed_reads} \\
        --threads ${params.threads} \\
        --very-sensitive \\
        | samtools view -bS -q 20 - \\
        | samtools sort -@ ${params.threads} -o ${sample_id}_aligned.bam -
    
    # Index BAM file
    samtools index ${sample_id}_aligned.bam
    
    # Generate alignment statistics
    samtools flagstat ${sample_id}_aligned.bam > ${sample_id}_alignment_stats.txt
    samtools stats ${sample_id}_aligned.bam >> ${sample_id}_alignment_stats.txt
    """
}

process PEAK_CALL {
    tag "${sample_id}"
    publishDir "${params.output_dir}/peaks", mode: 'copy'
    
    input:
    tuple val(sample_id), path(aligned_bam)
    
    output:
    tuple val(sample_id), path("${sample_id}_peaks.json"), emit: peaks_json
    tuple val(sample_id), path("${sample_id}_peaks.bed"), emit: peaks_bed
    tuple val(sample_id), path("${sample_id}_peak_stats.txt"), emit: peak_stats
    
    script:
    """
    # Call peaks using Rust implementation
    cargo run --release --bin peak_caller -- \\
        --input ${aligned_bam} \\
        --output ${sample_id}_peaks.json \\
        --window-size ${params.window_size} \\
        --min-coverage ${params.min_coverage} \\
        --pvalue-threshold ${params.pvalue_threshold} \\
        --threads ${params.threads} \\
        2>&1 | tee ${sample_id}_peak_calling.log
    
    # Convert JSON to BED format
    python3 -c "
import json
import sys

with open('${sample_id}_peaks.json', 'r') as f:
    peaks = json.load(f)

with open('${sample_id}_peaks.bed', 'w') as f:
    f.write('# chrom\\tstart\\tend\\tname\\tscore\\tstrand\\tsignalValue\\tpValue\\tqValue\\tpeak\\n')
    for i, peak in enumerate(peaks):
        score = min(1000, int(-10 * peak['pvalue'])) if peak['pvalue'] > 0 else 1000
        signal_value = peak['coverage']
        p_value = peak['pvalue']
        q_value = peak['qvalue']
        summit_offset = peak['summit'] - peak['start']
        
        f.write(f\"{peak['chrom']}\\t{peak['start']}\\t{peak['end']}\\t\"
                f\"peak_{i+1}\\t{score}\\t.\\t{signal_value:.2f}\\t\"
                f\"{p_value:.2e}\\t{q_value:.2e}\\t{summit_offset}\\n\")
    "
    
    # Generate peak statistics
    python3 -c "
import json
import statistics

with open('${sample_id}_peaks.json', 'r') as f:
    peaks = json.load(f)

if peaks:
    coverages = [p['coverage'] for p in peaks]
    pvalues = [p['pvalue'] for p in peaks]
    fold_enrichments = [p['fold_enrichment'] for p in peaks]
    
    stats = f'''Peak Calling Statistics for ${sample_id}
=================================
Total peaks: {len(peaks)}
Coverage statistics:
  Mean: {statistics.mean(coverages):.2f}
  Median: {statistics.median(coverages):.2f}
  Min: {min(coverages):.2f}
  Max: {max(coverages):.2f}
P-value statistics:
  Best (lowest): {min(pvalues):.2e}
  Worst (highest): {max(pvalues):.2e}
Fold enrichment statistics:
  Mean: {statistics.mean(fold_enrichments):.2f}
  Max: {max(fold_enrichments):.2f}
'''
    
    with open('${sample_id}_peak_stats.txt', 'w') as f:
        f.write(stats)
else:
    with open('${sample_id}_peak_stats.txt', 'w') as f:
        f.write('No peaks found for ${sample_id}\\n')
    "
    """
}

process MERGE_PEAKS {
    tag "merge_all_samples"
    publishDir "${params.output_dir}/merged", mode: 'copy'
    
    input:
    path peak_files
    
    output:
    path "merged_peaks.bed"
    path "merged_peak_summary.txt"
    
    script:
    """
    # Merge all BED files
    cat ${peak_files} | grep -v '^#' | sort -k1,1 -k2,2n > merged_peaks_raw.bed
    
    # Remove overlapping peaks (keep the one with highest score)
    bedtools merge -i merged_peaks_raw.bed -c 4,5,6,7,8,9,10 -o distinct,max,distinct,max,min,min,first > merged_peaks.bed
    
    # Generate summary statistics
    total_peaks=\$(wc -l < merged_peaks.bed)
    echo "Merged Peak Summary" > merged_peak_summary.txt
    echo "==================" >> merged_peak_summary.txt
    echo "Total merged peaks: \$total_peaks" >> merged_peak_summary.txt
    echo "" >> merged_peak_summary.txt
    echo "Peaks per chromosome:" >> merged_peak_summary.txt
    cut -f1 merged_peaks.bed | sort | uniq -c | sort -nr >> merged_peak_summary.txt
    """
}

process GENERATE_REPORT {
    tag "final_report"
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    path qc_reports
    path alignment_stats
    path peak_stats
    path merged_summary
    
    output:
    path "pipeline_report.html"
    
    script:
    """
    python3 -c "
import os
import glob
from datetime import datetime

html_content = '''<!DOCTYPE html>
<html>
<head>
    <title>Epigenomic Peak Calling Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        h1, h2 { color: #2c3e50; }
        table { border-collapse: collapse; width: 100%; margin: 20px 0; }
        th, td { border: 1px solid #ddd; padding: 12px; text-align: left; }
        th { background-color: #f2f2f2; }
        .summary { background-color: #e8f4f8; padding: 20px; border-radius: 5px; }
    </style>
</head>
<body>
    <h1>Epigenomic Peak Calling Pipeline Report</h1>
    <div class=\"summary\">
        <h2>Analysis Summary</h2>
        <p><strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        <p><strong>Pipeline Version:</strong> 1.0.0</p>
        <p><strong>Parameters:</strong></p>
        <ul>
            <li>Window Size: ${params.window_size} bp</li>
            <li>Min Coverage: ${params.min_coverage}</li>
            <li>P-value Threshold: ${params.pvalue_threshold}</li>
            <li>Threads: ${params.threads}</li>
        </ul>
    </div>
'''

# Add merged peak summary if available
try:
    with open('${merged_summary}', 'r') as f:
        summary_content = f.read()
    html_content += f'''
    <h2>Merged Peaks Summary</h2>
    <pre>{summary_content}</pre>
    '''
except:
    pass

# Add individual sample statistics
html_content += '<h2>Individual Sample Results</h2>'

for stats_file in glob.glob('*_peak_stats.txt'):
    sample_name = stats_file.replace('_peak_stats.txt', '')
    try:
        with open(stats_file, 'r') as f:
            stats_content = f.read()
        html_content += f'''
        <h3>{sample_name}</h3>
        <pre>{stats_content}</pre>
        '''
    except:
        continue

html_content += '''
</body>
</html>
'''

with open('pipeline_report.html', 'w') as f:
    f.write(html_content)
    "
    """
}

workflow {
    // Validate inputs first
    VALIDATE_INPUTS(params.input_dir, params.reference_genome)
    
    // Create input channel from FASTQ files
    fastq_ch = Channel
        .fromPath("${params.input_dir}/*.{fastq,fq}.gz")
        .map { file -> 
            def sample_id = file.baseName.toString().replaceAll(/\.(fastq|fq)$/, '')
            tuple(sample_id, file)
        }
    
    // Run FastQC for initial quality assessment
    QC_FASTQC(fastq_ch)
    
    // Quality control and trimming
    QC_AND_TRIM(fastq_ch)
    
    // Build reference index
    BUILD_INDEX(params.reference_genome)
    
    // Align reads
    ALIGN(QC_AND_TRIM.out.trimmed_reads, BUILD_INDEX.out.bt2_index)
    
    // Call peaks
    PEAK_CALL(ALIGN.out.aligned_bam)
    
    // Merge peaks from all samples
    MERGE_PEAKS(PEAK_CALL.out.peaks_bed.collect())
    
    // Generate final report
    GENERATE_REPORT(
        QC_AND_TRIM.out.qc_reports.collect(),
        ALIGN.out.stats.collect(),
        PEAK_CALL.out.peak_stats.collect(),
        MERGE_PEAKS.out[1]
    )
}

workflow.onComplete {
    log.info """
    ==============================================
    Pipeline completed at: ${new Date()}
    Success: ${workflow.success}
    Duration: ${workflow.duration}
    Work directory: ${workflow.workDir}
    Results directory: ${params.output_dir}
    ==============================================
    """.stripIndent()
}