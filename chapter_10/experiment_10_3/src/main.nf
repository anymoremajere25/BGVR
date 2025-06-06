nextflow.enable.dsl=2

params.bam_list       = 'bams.txt'
params.ref_intervals  = 'genome_intervals.txt'
params.parallel_chunck_size = 50000 // Example parameter for adjusting concurrency

// Read each line in bams.txt as a separate BAM file path
bamChannel = Channel.fromPath(params.bam_list)
    .splitText()
    .map { it.trim() }

// Read each line in genome_intervals.txt as a separate query interval (e.g., "chr1:1000-2000")
queriesChannel = Channel.fromPath(params.ref_intervals)
    .splitText()
    .map { it.trim() }

process coverageComputation {
    tag "${bam_file}"
    
    input:
    val(bam_file)

    output:
    path "coverage_${bam_file}.tsv"

    script:
    """
    rust_coverage_tool \
        --bam "${bam_file}" \
        --out "coverage_${bam_file}.tsv" \
        --chunk-size ${params.parallel_chunck_size}
    """
}

process mergeCoverage {
    input:
    path cov_files

    output:
    path "merged_coverage.tsv"

    script:
    """
    cat ${cov_files} > merged_coverage.tsv
    """
}

process intervalQuery {
    tag "${query_interval}"
    
    input:
    path merged_coverage
    val(query_interval)

    output:
    path "query_result_${query_interval}.tsv"

    script:
    """
    rust_interval_query_tool \
       --interval-file "${merged_coverage}" \
       --query "${query_interval}" \
       > "query_result_${query_interval}.tsv"
    """
}

workflow {
    // Compute coverage for each BAM
    coverage = coverageComputation(bamChannel)
    
    // Merge coverage files
    merged = mergeCoverage(coverage.collect())
    
    // Run interval queries in parallel
    intervalQuery(merged, queriesChannel)
}