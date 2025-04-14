#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.sample_list = 'samples.txt'

process collectQC {
    input:
    path bam_file

    output:
    path "qc_${bam_file.name}.json"

    """
    /home/ragon/BGVR/chapter_06/coverage_tool/target/release/coverage_tool \
      --bam ${bam_file} \
      --region chr1:1-300
    """
}

workflow {
    // Create channel from samples.txt
    def bamChannel = Channel.fromPath(params.sample_list).splitText().map { it.trim() }.map { file(it) }
    
    // Run collectQC process
    collectQC(bamChannel)
}
