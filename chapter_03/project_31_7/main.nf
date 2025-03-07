#!/usr/bin/env nextflow

/*
  A minimal Nextflow pipeline that:
    1) Compiles the Rust program (process 'compile')
    2) Executes the Rust binary on the provided FASTQ data (process 'analysis')

  Usage:
    nextflow run main.nf --fastq /path/to/reads.fastq
*/

params.fastq = params.fastq ?: 'example.fastq'
params.outdir = 'results'

process compile {
    /*
      This step compiles our Rust code so the next step can run it.
      In real pipelines, you might skip this if you already built
      your binary or if you store precompiled artifacts.
    */
    tag "Compiling code"
    // Ensure we have Rust & Cargo installed locally (e.g., via rustup)
    executor 'local'

    input:
    val dummy from Channel.of(1)  // a trivial input to trigger compilation

    output:
    file 'target/release/debruijn_bloom' into compiled_binary

    """
    echo "Compiling the Rust code..."
    cargo build --release
    cp target/release/debruijn_bloom .
    """
}

process analysis {
    /*
      This step runs the compiled Rust binary on a FASTQ file
      to build a de Bruijn graph and a Bloom filter.
    */
    tag "DeBruijn & Bloom"
    executor 'local'

    input:
    file 'debruijn_bloom' from compiled_binary
    path fastq_file from Channel.fromPath(params.fastq)

    output:
    file 'graph.json' into graph_out
    file 'bloom.json' into bloom_out

    """
    echo "Running the Rust analysis on ${fastq_file}..."
    ./debruijn_bloom --fastq ${fastq_file} --kmer 31 --outdir ${params.outdir}
    cp ${params.outdir}/graph.json ./
    cp ${params.outdir}/bloom.json ./
    """
}

workflow {
    // The pipeline simply compiles, then runs 'analysis' on the compiled binary.
    compile()
    analysis()
}
