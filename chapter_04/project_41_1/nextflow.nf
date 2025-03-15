#!/usr/bin/env nextflow

// Set a default value for the FASTQ file
params.synthetic_fastq = params.synthetic_fastq ?: 'synthetic_reads.fastq'

process buildPWMandMRF {
    input:
    path fastq_file from params.synthetic_fastq

    output:
    path 'pwm_results.txt'
    path 'mrf_results.txt'

    script:
    """
    cd rust_code
    cargo build --release
    ../rust_code/target/release/bioinformatics_tools ${fastq_file} pwm_results.txt mrf_results.txt
    """
}

workflow {
    buildPWMandMRF()
}
