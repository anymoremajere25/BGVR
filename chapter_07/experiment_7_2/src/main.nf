#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.ref_file      = 'reference.fa'
params.kmer_length   = 31
params.chunk_size    = 1000000
params.outdir        = 'results'
params.threads       = Runtime.runtime.availableProcessors()

// Change this to a mock command for testing - we'll simulate the tool functionality
params.kmer_tool     = "echo 'Mock kmer tool execution'"

process chunkReference {
    publishDir "${params.outdir}/chunks", mode: 'copy'

    input:
    path ref_path

    output:
    path "chunks/*.fa"

    script:
    """
    mkdir -p chunks
    # Use awk to split FASTA file into chunks
    awk '
    BEGIN {
        chunk_size=${params.chunk_size}
        overlap=${params.kmer_length - 1}
        chunk_num=1
        seq=""
        current_header=""
    }
    /^>/ {
        if (seq != "") {
            # Write previous sequence chunks
            seq_len = length(seq)
            for (i=1; i<=seq_len; i+=chunk_size) {
                end = i+chunk_size-1
                if (end > seq_len) end = seq_len
                
                # Add overlap unless this is first chunk
                start = i
                if (i > 1) start = i - overlap
                if (start < 1) start = 1
                
                chunk = substr(seq, start, end-start+1)
                print current_header " chunk:" chunk_num > ("chunks/chunk_" chunk_num ".fa")
                print chunk >> ("chunks/chunk_" chunk_num ".fa")
                chunk_num++
            }
        }
        current_header = \$0
        seq = ""
        next
    }
    {
        seq = seq \$0
    }
    END {
        # Process the last sequence
        if (seq != "") {
            seq_len = length(seq)
            for (i=1; i<=seq_len; i+=chunk_size) {
                end = i+chunk_size-1
                if (end > seq_len) end = seq_len
                
                # Add overlap unless this is first chunk
                start = i
                if (i > 1) start = i - overlap
                if (start < 1) start = 1
                
                chunk = substr(seq, start, end-start+1)
                print current_header " chunk:" chunk_num > ("chunks/chunk_" chunk_num ".fa")
                print chunk >> ("chunks/chunk_" chunk_num ".fa")
                chunk_num++
            }
        }
    }' ${ref_path}
    """
}

process buildPartialIndex {
    cpus params.threads
    publishDir "${params.outdir}/partial", mode: 'copy'

    input:
    path chunk

    output:
    path "partial_*.json"

    script:
    def outname = "partial_${chunk.baseName}.json"
    """
    # Create a mock JSON output for testing since the real tool isn't accessible
    echo '{
      "kmers": {
        "ACGT": 5,
        "CGTG": 3,
        "GTAC": 2
      },
      "metadata": {
        "filename": "${chunk}",
        "k": ${params.kmer_length},
        "total_kmers": 10
      }
    }' > ${outname}
    
    # Log what would have been executed
    echo "Would execute: ${params.kmer_tool} --input ${chunk} --output ${outname}" > ${outname}.log
    """
}

process mergeIndexes {
    publishDir params.outdir, mode: 'copy'

    input:
    path partial_files

    output:
    path "global_index.json"

    script:
    def all_files = partial_files.join(" ")
    """
    # Create a mock merged JSON output for testing
    echo '{
      "kmers": {
        "ACGT": 10,
        "CGTG": 6,
        "GTAC": 4,
        "TACG": 3
      },
      "metadata": {
        "input_files": ${partial_files.size()},
        "k": ${params.kmer_length},
        "total_kmers": 23
      }
    }' > global_index.json
    
    # Log what would have been executed
    echo "Would execute: ${params.kmer_tool} --merge ${all_files} --output global_index.json" > merge_command.log
    """
}

workflow {
    // Create a channel for the input reference file
    ref_file_channel = Channel.fromPath(params.ref_file)

    // Run chunkReference process
    chunk_channel = chunkReference(ref_file_channel)

    // Run buildPartialIndex process with chunked files as input
    index_channel = buildPartialIndex(chunk_channel.flatten())

    // Run mergeIndexes process using the output from buildPartialIndex
    mergeIndexes(index_channel.collect())
}