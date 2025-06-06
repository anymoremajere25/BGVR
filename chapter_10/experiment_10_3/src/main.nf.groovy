#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Comprehensive QC and Signal Normalization Pipeline Parameters
params.input_dir = './data'
params.output_dir = './results'
params.reference_genome = './reference/genome.fa'
params.background_regions = null  // Optional background regions BED
params.blacklist_regions = null   // Optional blacklist regions BED
params.control_samples = null     // Optional control/input samples directory

// Quality control parameters
params.min_mapq = 10
params.fragment_shift = 75
params.extend_reads = 150
params.bin_size = 1000
params.peak_window_size = 1000
params.min_peak_coverage = 5.0

// Normalization parameters
params.normalize_rpm = true
params.normalize_rpkm = true
params.normalize_tpm = false
params.generate_bigwig = true

// Analysis parameters
params.threads = 4
params.memory = '8 GB'
params.time_limit = '4h'

// Quality thresholds
params.min_frip = 0.01
params.min_snr = 2.0
params.min_complexity = 0.7
params.max_duplication = 0.5

log.info """
===============================================
Signal Normalization and Quality Control Pipeline
===============================================
Input directory      : ${params.input_dir}
Output directory     : ${params.output_dir}
Reference genome     : ${params.reference_genome}
Background regions   : ${params.background_regions ?: 'Genome-wide estimate'}
Blacklist regions    : ${params.blacklist_regions ?: 'None'}
Control samples      : ${params.control_samples ?: 'None'}

QC Parameters:
- Min mapping quality: ${params.min_mapq}
- Fragment shift     : ${params.fragment_shift} bp
- Read extension     : ${params.extend_reads} bp
- Bin size          : ${params.bin_size} bp

Quality Thresholds:
- Min FRiP score    : ${params.min_frip}
- Min SNR           : ${params.min_snr}
- Min complexity    : ${params.min_complexity}
- Max duplication   : ${params.max_duplication}

Threads: ${params.threads}
===============================================
"""

// Validate inputs
process VALIDATE_QC_INPUTS {
    tag "qc_validation"
    
    input:
    path input_dir
    path reference_genome
    path background_regions
    path blacklist_regions
    
    output:
    stdout
    
    script:
    def bg_check = background_regions.name != 'NO_FILE' ? 
        "test -f ${background_regions} || (echo 'ERROR: Background regions file not found' && exit 1)" :
        "echo 'Using genome-wide background estimation'"
    
    def bl_check = blacklist_regions.name != 'NO_FILE' ? 
        "test -f ${blacklist_regions} || (echo 'ERROR: Blacklist regions file not found' && exit 1)" :
        "echo 'No blacklist regions specified'"
    
    """
    echo "=== QC Pipeline Input Validation ==="
    
    # Check input directory
    test -d "${input_dir}" || (echo "ERROR: Input directory not found: ${input_dir}" && exit 1)
    
    # Check reference genome
    test -f "${reference_genome}" || (echo "ERROR: Reference genome not found: ${reference_genome}" && exit 1)
    
    # Check optional files
    ${bg_check}
    ${bl_check}
    
    # Count FASTQ files
    fastq_count=\$(find ${input_dir} -name "*.fastq.gz" -o -name "*.fq.gz" | wc -l)
    test \$fastq_count -gt 0 || (echo "ERROR: No FASTQ files found in ${input_dir}" && exit 1)
    
    echo "✓ Found \$fastq_count FASTQ files"
    echo "✓ All inputs validated successfully"
    """
}

// Build comprehensive reference indices
process BUILD_QC_INDICES {
    tag "qc_indexing"
    publishDir "${params.output_dir}/reference", mode: 'copy'
    
    input:
    path reference_genome
    
    output:
    path "genome_index*", emit: bt2_index
    path "${reference_genome}.fai", emit: fai_index
    path "chrom_sizes.txt", emit: chrom_sizes
    
    script:
    """
    echo "Building indices for QC pipeline..."
    
    # Build bowtie2 index
    bowtie2-build --threads ${params.threads} ${reference_genome} genome_index
    
    # Build samtools index
    samtools faidx ${reference_genome}
    
    # Create chromosome sizes file for normalization
    cut -f1,2 ${reference_genome}.fai > chrom_sizes.txt
    
    echo "Reference indexing completed"
    """
}

// Enhanced read trimming and QC
process ENHANCED_QC_AND_TRIM {
    tag "${sample_id}"
    publishDir "${params.output_dir}/qc/trimmed", mode: 'copy', pattern: "*.{json,html}"
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq.gz"), emit: trimmed_reads
    tuple val(sample_id), path("${sample_id}_fastp.{json,html}"), emit: qc_reports
    tuple val(sample_id), path("${sample_id}_pre_alignment_stats.txt"), emit: pre_stats
    
    script:
    """
    # Enhanced QC and trimming
    fastp \\
        -i ${reads} \\
        -o ${sample_id}_trimmed.fastq.gz \\
        --qualified_quality_phred 20 \\
        --unqualified_percent_limit 20 \\
        --length_required 36 \\
        --low_complexity_filter \\
        --complexity_threshold 30 \\
        --detect_adapter_for_pe \\
        --thread ${params.threads} \\
        --json ${sample_id}_fastp.json \\
        --html ${sample_id}_fastp.html
    
    # Extract pre-alignment statistics
    echo "Sample: ${sample_id}" > ${sample_id}_pre_alignment_stats.txt
    echo "Input reads: \$(jq '.summary.before_filtering.total_reads // 0' ${sample_id}_fastp.json)" >> ${sample_id}_pre_alignment_stats.txt
    echo "Output reads: \$(jq '.summary.after_filtering.total_reads // 0' ${sample_id}_fastp.json)" >> ${sample_id}_pre_alignment_stats.txt
    echo "Read length: \$(jq '.summary.after_filtering.read1_mean_length // 0' ${sample_id}_fastp.json)" >> ${sample_id}_pre_alignment_stats.txt
    echo "GC content: \$(jq '.summary.after_filtering.gc_content // 0' ${sample_id}_fastp.json)" >> ${sample_id}_pre_alignment_stats.txt
    echo "Q30 rate: \$(jq '.summary.after_filtering.q30_rate // 0' ${sample_id}_fastp.json)" >> ${sample_id}_pre_alignment_stats.txt
    
    # Calculate adapter content
    adapter_rate=\$(jq '.adapter_cutting.adapter_trimmed_reads // 0' ${sample_id}_fastp.json)
    total_reads=\$(jq '.summary.before_filtering.total_reads // 1' ${sample_id}_fastp.json)
    adapter_percent=\$(echo "scale=3; \$adapter_rate * 100 / \$total_reads" | bc -l)
    echo "Adapter content: \$adapter_percent%" >> ${sample_id}_pre_alignment_stats.txt
    """
}

// Enhanced alignment with QC-specific parameters
process QC_ALIGNMENT {
    tag "${sample_id}"
    publishDir "${params.output_dir}/alignments", mode: 'copy', pattern: "*.{bam,bai}"
    
    input:
    tuple val(sample_id), path(trimmed_reads)
    path bt2_index
    
    output:
    tuple val(sample_id), path("${sample_id}_aligned.bam"), emit: aligned_bam
    tuple val(sample_id), path("${sample_id}_aligned.bam.bai"), emit: bam_index
    tuple val(sample_id), path("${sample_id}_alignment_detailed.txt"), emit: alignment_stats
    
    script:
    """
    echo "Aligning ${sample_id} with QC-optimized parameters..."
    
    # Alignment optimized for QC analysis
    bowtie2 \\
        -x genome_index \\
        -U ${trimmed_reads} \\
        --threads ${params.threads} \\
        --very-sensitive \\
        --no-discordant \\
        --no-mixed \\
        -k 1 \\
        2> ${sample_id}_bowtie2.log \\
        | samtools view -bS -F 4 - \\
        | samtools sort -@ ${params.threads} -o ${sample_id}_aligned.bam -
    
    # Index BAM file
    samtools index ${sample_id}_aligned.bam
    
    # Generate detailed alignment statistics for QC
    samtools flagstat ${sample_id}_aligned.bam > ${sample_id}_alignment_detailed.txt
    samtools stats ${sample_id}_aligned.bam >> ${sample_id}_alignment_detailed.txt
    
    # Extract key metrics
    echo "" >> ${sample_id}_alignment_detailed.txt
    echo "=== Key QC Metrics ===" >> ${sample_id}_alignment_detailed.txt
    echo "Total reads: \$(samtools view -c ${sample_id}_aligned.bam)" >> ${sample_id}_alignment_detailed.txt
    echo "Mapped reads: \$(samtools view -c -F 4 ${sample_id}_aligned.bam)" >> ${sample_id}_alignment_detailed.txt
    echo "Properly paired: \$(samtools view -c -f 2 ${sample_id}_aligned.bam)" >> ${sample_id}_alignment_detailed.txt
    echo "Duplicates: \$(samtools view -c -f 1024 ${sample_id}_aligned.bam)" >> ${sample_id}_alignment_detailed.txt
    echo "High quality (MAPQ>=30): \$(samtools view -c -q 30 ${sample_id}_aligned.bam)" >> ${sample_id}_alignment_detailed.txt
    
    # Calculate mapping rate
    total=\$(samtools view -c ${sample_id}_aligned.bam)
    mapped=\$(samtools view -c -F 4 ${sample_id}_aligned.bam)
    if [ \$total -gt 0 ]; then
        mapping_rate=\$(echo "scale=3; \$mapped * 100 / \$total" | bc -l)
        echo "Mapping rate: \$mapping_rate%" >> ${sample_id}_alignment_detailed.txt
    fi
    """
}

// Simple peak calling for QC analysis
process SIMPLE_PEAK_CALLING {
    tag "${sample_id}"
    publishDir "${params.output_dir}/peaks", mode: 'copy'
    
    input:
    tuple val(sample_id), path(aligned_bam)
    
    output:
    tuple val(sample_id), path("${sample_id}_peaks.bed"), emit: peaks_bed
    tuple val(sample_id), path("${sample_id}_peak_calling.log"), emit: peak_log
    
    script:
    """
    echo "Calling peaks for QC analysis: ${sample_id}"
    
    # Use Rust peak caller
    cargo run --release --bin rust_peak_caller -- \\
        ${aligned_bam} \\
        --output ${sample_id}_peaks.bed \\
        --window-size ${params.peak_window_size} \\
        --min-coverage ${params.min_peak_coverage} \\
        --threads ${params.threads} \\
        2>&1 | tee ${sample_id}_peak_calling.log
    
    # Ensure output file exists even if no peaks found
    if [ ! -f "${sample_id}_peaks.bed" ]; then
        echo "# No peaks found" > ${sample_id}_peaks.bed
    fi
    
    # Add peak count to log
    peak_count=\$(grep -v "^#" ${sample_id}_peaks.bed | wc -l)
    echo "Total peaks found: \$peak_count" >> ${sample_id}_peak_calling.log
    """
}

// Generate background regions if not provided
process GENERATE_BACKGROUND_REGIONS {
    tag "background_generation"
    publishDir "${params.output_dir}/regions", mode: 'copy'
    
    input:
    path chrom_sizes
    
    output:
    path "background_regions.bed", emit: background_bed
    
    when:
    !params.background_regions
    
    script:
    """
    echo "Generating random background regions..."
    
    # Generate random background regions (1% of genome in 1kb windows)
    python3 << 'EOF'
import random
import sys

# Read chromosome sizes
chrom_sizes = {}
with open('${chrom_sizes}', 'r') as f:
    for line in f:
        chrom, size = line.strip().split()
        chrom_sizes[chrom] = int(size)

# Generate background regions
window_size = ${params.bin_size}
total_genome_size = sum(chrom_sizes.values())
num_regions = max(1000, total_genome_size // (window_size * 100))  # ~1% coverage

background_regions = []
for _ in range(num_regions):
    chrom = random.choice(list(chrom_sizes.keys()))
    chrom_size = chrom_sizes[chrom]
    
    if chrom_size > window_size * 2:
        start = random.randint(0, chrom_size - window_size)
        end = start + window_size
        background_regions.append((chrom, start, end))

# Write to BED file
with open('background_regions.bed', 'w') as f:
    f.write("# Randomly generated background regions\\n")
    for i, (chrom, start, end) in enumerate(background_regions):
        f.write(f"{chrom}\\t{start}\\t{end}\\tbackground_{i+1}\\n")

print(f"Generated {len(background_regions)} background regions")
EOF
    """
}

// Comprehensive quality control analysis
process COMPREHENSIVE_QC_ANALYSIS {
    tag "${sample_id}"
    publishDir "${params.output_dir}/qc/metrics", mode: 'copy'
    
    input:
    tuple val(sample_id), path(aligned_bam), path(peaks_bed)
    path background_regions
    path blacklist_regions
    
    output:
    tuple val(sample_id), path("${sample_id}_qc_results.json"), emit: qc_json
    tuple val(sample_id), path("${sample_id}_qc_summary.txt"), emit: qc_summary
    tuple val(sample_id), path("${sample_id}_qc_detailed.log"), emit: qc_log
    
    script:
    def bg_option = background_regions.name != 'NO_FILE' ? "--background-regions ${background_regions}" : ""
    def bl_option = blacklist_regions.name != 'NO_FILE' ? "--blacklist-regions ${blacklist_regions}" : ""
    
    """
    echo "Running comprehensive QC analysis for ${sample_id}..."
    
    # Run Rust QC analysis
    cargo run --release --bin rust_qc -- \\
        --input ${aligned_bam} \\
        --peaks ${peaks_bed} \\
        --output ${sample_id}_qc_results.json \\
        ${bg_option} \\
        ${bl_option} \\
        --fragment-shift ${params.fragment_shift} \\
        --extend-reads ${params.extend_reads} \\
        --bin-size ${params.bin_size} \\
        --min-mapq ${params.min_mapq} \\
        --threads ${params.threads} \\
        --verbose \\
        2>&1 | tee ${sample_id}_qc_detailed.log
    
    # Generate human-readable summary
    python3 << 'EOF'
import json
import sys

try:
    with open('${sample_id}_qc_results.json', 'r') as f:
        qc_data = json.load(f)
    
    # Create summary report
    summary = f"""Quality Control Summary for ${sample_id}
{'=' * 50}
Date: $(date)

Basic Statistics:
- Total reads: {qc_data.get('total_reads', 'N/A'):,}
- Mapped reads: {qc_data.get('mapped_reads', 'N/A'):,}
- Mapping rate: {qc_data.get('mapping_rate', 0)*100:.1f}%
- Duplication rate: {qc_data.get('duplication_rate', 0)*100:.1f}%
- Library complexity: {qc_data.get('library_complexity', 0):.3f}

Signal Quality:
- FRiP score: {qc_data.get('frip_score', 0):.3f}
- Signal-to-noise ratio: {qc_data.get('signal_to_noise_ratio', 0):.2f}
- Reads in peaks: {qc_data.get('reads_in_peaks', 'N/A'):,}

Fragment Analysis:
- Mean fragment size: {qc_data.get('mean_fragment_size', 0):.1f} bp
- Fragment size std: {qc_data.get('fragment_size_std', 0):.1f} bp

Normalization Factors:
- RPM factor: {qc_data.get('rpm_factor', 0):.6f}
- RPKM factor: {qc_data.get('rpkm_factor', 0):.6f}

Quality Thresholds:
- FRiP threshold (>={params.min_frip}): {'PASS' if qc_data.get('pass_frip_threshold', False) else 'FAIL'}
- SNR threshold (>={params.min_snr}): {'PASS' if qc_data.get('pass_snr_threshold', False) else 'FAIL'}
- Complexity threshold (>={params.min_complexity}): {'PASS' if qc_data.get('pass_complexity_threshold', False) else 'FAIL'}

Overall Quality: {qc_data.get('overall_quality', 'UNKNOWN')}
"""
    
    with open('${sample_id}_qc_summary.txt', 'w') as f:
        f.write(summary)
    
    print("QC summary generated successfully")

except Exception as e:
    print(f"Error generating QC summary: {e}")
    with open('${sample_id}_qc_summary.txt', 'w') as f:
        f.write(f"Error generating QC summary for ${sample_id}: {e}\\n")
EOF
    """
}

// Signal normalization and BigWig generation
process SIGNAL_NORMALIZATION {
    tag "${sample_id}"
    publishDir "${params.output_dir}/normalized", mode: 'copy'
    
    input:
    tuple val(sample_id), path(aligned_bam), path(qc_json)
    path chrom_sizes
    
    output:
    tuple val(sample_id), path("${sample_id}_rpm.bw"), emit: rpm_bigwig, optional: true
    tuple val(sample_id), path("${sample_id}_rpkm.bw"), emit: rpkm_bigwig, optional: true
    tuple val(sample_id), path("${sample_id}_normalized.bedgraph"), emit: bedgraph
    
    when:
    params.normalize_rpm || params.normalize_rpkm || params.generate_bigwig
    
    script:
    """
    echo "Generating normalized signal tracks for ${sample_id}..."
    
    # Extract normalization factor from QC results
    rpm_factor=\$(jq '.rpm_factor // 1.0' ${qc_json})
    
    # Generate coverage bedgraph
    bedtools genomecov \\
        -ibam ${aligned_bam} \\
        -bg \\
        -scale \$rpm_factor \\
        > ${sample_id}_rpm.bedgraph
    
    # Sort bedgraph
    sort -k1,1 -k2,2n ${sample_id}_rpm.bedgraph > ${sample_id}_normalized.bedgraph
    
    # Generate BigWig files if requested
    if [ "${params.generate_bigwig}" = "true" ]; then
        # Convert to BigWig
        bedGraphToBigWig ${sample_id}_normalized.bedgraph ${chrom_sizes} ${sample_id}_rpm.bw
        
        # Generate RPKM normalized version (placeholder - would need gene annotations for proper RPKM)
        cp ${sample_id}_rpm.bw ${sample_id}_rpkm.bw
    fi
    
    echo "Signal normalization completed for ${sample_id}"
    """
}

// Aggregate QC results across all samples
process AGGREGATE_QC_RESULTS {
    tag "qc_aggregation"
    publishDir "${params.output_dir}/summary", mode: 'copy'
    
    input:
    path qc_jsons
    path qc_summaries
    
    output:
    path "qc_aggregate_report.html", emit: aggregate_html
    path "qc_metrics_table.csv", emit: metrics_csv
    path "qc_pass_fail_summary.txt", emit: pass_fail_summary
    
    script:
    """
    echo "Aggregating QC results across all samples..."
    
    python3 << 'EOF'
import json
import pandas as pd
import os
from datetime import datetime

# Collect all QC data
qc_data = []
for qc_file in [f for f in os.listdir('.') if f.endswith('_qc_results.json')]:
    sample_id = qc_file.replace('_qc_results.json', '')
    try:
        with open(qc_file, 'r') as f:
            data = json.load(f)
            data['sample_id'] = sample_id
            qc_data.append(data)
    except Exception as e:
        print(f"Error loading {qc_file}: {e}")

if not qc_data:
    print("No QC data found")
    # Create empty outputs
    with open('qc_aggregate_report.html', 'w') as f:
        f.write('<html><body><h1>No QC data available</h1></body></html>')
    with open('qc_metrics_table.csv', 'w') as f:
        f.write('sample_id,error\\n')
    with open('qc_pass_fail_summary.txt', 'w') as f:
        f.write('No samples processed\\n')
    exit(0)

# Create DataFrame
df = pd.DataFrame(qc_data)

# Generate CSV table
csv_columns = [
    'sample_id', 'total_reads', 'mapping_rate', 'duplication_rate',
    'frip_score', 'signal_to_noise_ratio', 'library_complexity',
    'mean_fragment_size', '