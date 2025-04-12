params.bam_list = 'bams.txt'
params.ref_intervals = 'genome_intervals.txt'
params.parallel_chunk_size = 50000

process coverageComputation {
    tag "${bam_file}"
    input:
    path bam_file
    output:
    path "coverage_${bam_file}.tsv"
    """
    /home/ragon/BGVR/chapter_06/experiment_6_5_1_b/target/debug/genomic_interval_tree coverage \
        --bam ${bam_file} \
        --out coverage_${bam_file}.tsv \
        --chunk-size ${params.parallel_chunk_size}
    """
}

process mergeCoverage {
    input:
    path cov_files
    output:
    path "merged_coverage.tsv"
    """
    cat ${cov_files} > merged_coverage.tsv
    """
}

process intervalQuery {
    tag "${query_interval}"
    input:
    path intervals
    val query_interval
    output:
    path "query_result_${query_interval}.tsv"
    """
    /home/ragon/BGVR/chapter_06/experiment_6_5_1_b/target/debug/genomic_interval_tree query \
        --interval-file ${intervals} \
        --query ${query_interval} \
        --output query_result_${query_interval}.tsv
    """
}

workflow {
    bams = Channel.fromPath(params.bam_list).splitText().map { it.trim() }.map { file(it) }
    queries = Channel.fromPath(params.ref_intervals).splitText().map { it.trim() }
    coverage_files = coverageComputation(bams)
    merged_intervals = mergeCoverage(coverage_files.collect())
    intervalQuery(merged_intervals, queries)
}