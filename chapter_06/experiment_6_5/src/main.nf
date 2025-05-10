nextflow.enable.dsl=2

params.bam_list       = 'bams.txt'
params.ref_intervals  = 'genome_intervals.txt'
params.parallel_chunck_size = 50000 // Example parameter for adjusting concurrency
params.mock = true // Set to true to use mock commands for testing

workflow {
    // Define the channels for BAM files and reference intervals
    bamChannel = Channel.fromPath(params.bam_list).splitText().map { it.trim() }
    queriesChannel = Channel.fromPath(params.ref_intervals).splitText().map { it.trim() }
    
    // Run the processes in the workflow
    coverageResults = coverageComputation(bamChannel)
    mergedCoverage = mergeCoverage(coverageResults.collect())
    intervalQuery(mergedCoverage, queriesChannel)
}

process coverageComputation {
    tag "${bam_file}"
    
    input:
    val bam_file
    
    output:
    path "coverage_${bam_file}.tsv"
    
    script:
    if (params.mock) {
        """
        # Mock coverage computation
        echo "Creating mock coverage file for ${bam_file}"
        echo "chrom\tstart\tend\tdepth\tsample" > coverage_${bam_file}.tsv
        echo "chr1\t100\t200\t10\t${bam_file}" >> coverage_${bam_file}.tsv
        echo "chr1\t200\t300\t15\t${bam_file}" >> coverage_${bam_file}.tsv
        echo "chr1\t300\t400\t20\t${bam_file}" >> coverage_${bam_file}.tsv
        echo "chr2\t100\t200\t5\t${bam_file}" >> coverage_${bam_file}.tsv
        echo "chr2\t200\t300\t8\t${bam_file}" >> coverage_${bam_file}.tsv
        """
    } else {
        """
        # The Rust coverage tool here would parse the BAM file via rust-htslib,
        # compute coverage intervals (possibly using an ndarray for efficient numeric ops),
        # and output them as TSV. Each ephemeral task handles one BAM.
        /mnt/c/Users/trian/BGVR/chapter_06/experiment_65/target/debug/rust_coverage_tool \
            --bam ${bam_file} \
            --out coverage_${bam_file}.tsv \
            --chunk-size ${params.parallel_chunck_size}
        """
    }
}

process mergeCoverage {
    input:
    path cov_files
    
    output:
    path "merged_coverage.tsv"
    
    script:
    if (params.mock) {
        """
        # Mock merge operation
        echo "Merging coverage files: ${cov_files}"
        # Get header from first file
        head -n 1 \$(ls ${cov_files} | head -n 1) > merged_coverage.tsv
        # Append data from all files (skipping headers)
        for file in ${cov_files}; do
            tail -n +2 \$file >> merged_coverage.tsv
        done
        """
    } else {
        """
        # Merge individual coverage files into one TSV, ready for interval-based queries.
        # Real-world pipelines often sort or index these data for faster lookups.
        cat ${cov_files} > merged_coverage.tsv
        """
    }
}

process intervalQuery {
    tag "${query_interval}"
    
    input:
    path merged_coverage
    val query_interval
    
    output:
    path "query_result_${query_interval}.tsv"
    
    script:
    if (params.mock) {
        """
        # Mock interval query operation
        echo "Querying interval: ${query_interval}"
        # Parse query interval
        QUERY=\$(echo "${query_interval}" | tr -d '\\n')
        IFS=':' read -r CHROM RANGE <<< "\$QUERY"
        IFS='-' read -r START END <<< "\$RANGE"
        
        # Create header
        echo "chrom\tstart\tend\tdepth\tsample" > query_result_${query_interval}.tsv
        
        # Simple grep to simulate interval overlap (in real tool, this would use proper interval tree)
        grep "\$CHROM" ${merged_coverage} | while read -r line; do
            echo "\$line" >> query_result_${query_interval}.tsv
        done
        """
    } else {
        """
        # The Rust interval query tool reads merged coverage intervals, builds an interval tree,
        # and retrieves overlapping intervals for the specified query range. This can occur in parallel
        # when multiple queries are submitted concurrently in HPC or cloud environments.
        /mnt/c/Users/trian/BGVR/chapter_06/experiment_65/target/debug/interval_query_tool \
           --interval-file ${merged_coverage} \
           --query "${query_interval}" \
           > query_result_${query_interval}.tsv
        """
    }
}