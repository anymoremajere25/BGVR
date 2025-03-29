#!/usr/bin/env nextflow

/*
 * main.nf
 *
 * Example usage:
 *    nextflow run main.nf --input_fasta genome.fa --chunk_size 1000000
 * 
 * This pipeline:
 *  1) Splits the input genome FASTA into smaller chunks.
 *  2) Runs the Rust motif-scanning code on each chunk (in parallel).
 *  3) Merges the partial JSON outputs into one consolidated file.
 */

params.input_fasta   = "genome.fa"
params.chunk_size    = 100000  // chunk size in bp
params.output_prefix = "tata_scan"

// Define the Rust project directory (parent of current directory)
// Fix: Use toAbsolutePath() for path resolution
rustProjectDir = projectDir.parent.toAbsolutePath().toString()

/*
 * Process 1: Chunk the large FASTA into smaller files for parallel processing.
 */
process splitFasta {
    executor = 'local'
    
    input:
    path 'input.fa'
    
    output:
    path '*.fa', emit: chunks
    
    script:
    """
    # For large HPC usage, you might prefer a more advanced tool like 'pyfaidx' or 'seqkit' to chunk by size.
    seqkit split2 --by-size ${params.chunk_size} --out-dir . input.fa
    # This command typically produces multiple files named input.part_001.fa, input.part_002.fa, etc.
    """
}

/*
 * Process 2: Run the Rust code on each chunk, producing partial JSON hits.
 */
process scanMotif {
    executor = 'local'
    
    input:
    path chunk
    
    output:
    path "partial_*.json", emit: partial_jsons
    
    script:
    """
    # Get current working directory
    WORK_DIR=\$PWD
    
    # Navigate to the project directory where Cargo.toml is located
    cd ${rustProjectDir}
    
    # Build the Rust project in the parent directory
    cargo build --release
    
    # Debug output
    echo "Rust project directory: ${rustProjectDir}"
    echo "Working directory: \$WORK_DIR"
    echo "Chunk file: \$WORK_DIR/${chunk}"
    
    # Run the binary with the full path to the chunk file
    ${rustProjectDir}/target/release/experiment_43_5 \$WORK_DIR/${chunk} \$WORK_DIR/partial_${chunk.simpleName}.json
    """
}

/*
 * Process 3: Merge all partial JSON results into one file.
 */
process mergeHits {
    executor = 'local'
    publishDir '.', mode: 'copy'
    
    input:
    path partials
    
    output:
    path "${params.output_prefix}_merged.json"
    
    script:
    """
    # A simple merge by concatenating JSON lines
    cat ${partials} > ${params.output_prefix}_merged.json
    """
}

workflow {
    input_ch = channel.fromPath(params.input_fasta)
    
    splitFasta(input_ch)
    scanMotif(splitFasta.out.chunks.flatten())
    mergeHits(scanMotif.out.partial_jsons.collect())
}