#!/usr/bin/env nextflow

/*
  A simple Nextflow pipeline:
    1) Compiles the Rust code (process 'compile').
    2) Executes the Rust binary on the provided FASTA file (process 'analysis').

  Usage:
    nextflow run main.nf --fasta /path/to/large_sequence.fa
*/

params.fasta = params.fasta ?: 'large_sequence.fa'
params.outfile = 'partial_suffix_arrays.json'

process compile {
    executor 'local'

    input:
    val dummy from Channel.of(1)  // Ensures the process runs once

    output:
    path 'partial_sa' into compiled_binary  // Use `path` instead of `file`

    """
    echo "Compiling the Rust code..."
    cargo build --release
    cp target/release/project_32_3 partial_sa
    """
}

process analysis {
    executor 'local'

    input:
    path partial_sa from compiled_binary
    path fasta_file from params.fasta

    output:
    path params.outfile into results

    """
    echo "Running partial suffix array on ${fasta_file}..."
    ./partial_sa ${fasta_file}
    cp partial_suffix_arrays.json ${params.outfile}
    """
}

workflow {
    compile()
    analysis()
}

