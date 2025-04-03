nextflow.enable.dsl=2

params.bam_list       = 'bams.txt'
params.ref_intervals  = 'genome_intervals.txt'
params.parallel_chunck_size = 50000 // Example parameter for adjusting concurrency

// Read each line in bams.txt as a separate BAM file path
bamChannel = Channel.fromPath(params.bam_list).splitText()

// Read each line in genome_intervals.txt as a separate query interval (e.g., "chr1:1000-2000")
queriesChannel = Channel.fromPath(params.ref_intervals).splitText()

process coverageComputation {
    tag "${bam_file}"
    
    input:
    val bam_file from bamChannel

    output:
    file "coverage_${bam_file}.tsv"

    script:
    """
    # The Rust coverage tool here would parse the BAM file via rust-htslib,
    # compute coverage intervals (possibly using an ndarray for efficient numeric ops),
    # and output them as TSV. Each ephemeral task handles one BAM.
    rust_coverage_tool \
        --bam ${bam_file} \
        --out coverage_${bam_file}.tsv \
        --chunk-size ${params.parallel_chunck_size}
    """
}

process mergeCoverage {
    input:
    file cov_files from coverageComputation.out.collect()

    output:
    file "merged_coverage.tsv"

    script:
    """
    # Merge individual coverage files into one TSV, ready for interval-based queries.
    # Real-world pipelines often sort or index these data for faster lookups.
    cat coverage_*.tsv > merged_coverage.tsv
    """
}

process intervalQuery {
    tag "${query_interval}"
    
    input:
    file "merged_coverage.tsv"
    val query_interval from queriesChannel

    output:
    file "query_result_${query_interval}.tsv"

    script:
    """
    # The Rust interval query tool reads merged coverage intervals, builds an interval tree,
    # and retrieves overlapping intervals for the specified query range. This can occur in parallel
    # when multiple queries are submitted concurrently in HPC or cloud environments.
    rust_interval_query_tool \
       --interval-file merged_coverage.tsv \
       --query ${query_interval} \
       > query_result_${query_interval}.tsv
    """
}

workflow {
    // Compute coverage for each BAM, merge results, and run interval queries in parallel
    coverageComputation(bamChannel)
    mergeCoverage(coverageComputation.out)
    intervalQuery(mergeCoverage.out, queriesChannel)
}
