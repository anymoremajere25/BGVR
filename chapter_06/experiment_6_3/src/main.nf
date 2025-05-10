params.input_bcf   = '/mnt/c/Users/trian/BGVR/chapter_06/experiment_63/src/wgs_cohort.bcf'
params.chunks_file = 'chunks_list.txt'
params.min_qual    = 30
params.min_depth   = 10
params.chunk_size  = 50000  // optional tuning

process filterVariants {
    input:
    path bcf_file
    val chunk

    output:
    file("filtered_${chunk.replace(':','_')}.bcf")

    script:
    def out_file = "filtered_${chunk.replace(':','_')}.bcf"
    def bcf_filter_tool_path = '/mnt/c/Users/trian/BGVR/chapter_06/experiment_63/target/debug/bcf_filter_tool'

    """
    ${bcf_filter_tool_path} \\
        --input ${bcf_file} \\
        --output ${out_file} \\
        --min-qual ${params.min_qual} \\
        --min-depth ${params.min_depth} \\
        --chunk-size ${params.chunk_size}
    """
}

workflow {
    // Defining input channels
    Channel
        .fromPath(params.input_bcf)
        .set { bcf_file_channel }

    Channel
        .from(['chr1:1000-1200', 'chr2:2000-2200', 'chr3:3000-3200'])
        .set { chunk_channel }

    // Passing the combined channels to the filterVariants process directly within the workflow context
    filterVariants(bcf_file_channel, chunk_channel)
        .view()  // Optionally visualize the output
}
