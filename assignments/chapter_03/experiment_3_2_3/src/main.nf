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
    // This step compiles our Rust code using Cargo
    executor 'local'

    input:
    val dummy from Channel.of(1)  // a trivial input to trigger compilation

    output:
    file 'target/release/partial_sa' into compiled_binary

    """
    echo "Compiling the Rust code..."
    cargo build --release
    cp target/release/partial_sa .
    """
}

process analysis {
    // This step runs the Rust binary on the FASTA file to build partial suffix arrays
    executor 'local'

    input:
    file 'partial_sa' from compiled_binary
    path fasta_file from Channel.fromPath(params.fasta)

    output:
    file(params.outfile) into results

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
