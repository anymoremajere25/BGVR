## 3.1. Introduction to Data Processing in Bioinformatics

## project_31_7
This Rust and Nextflow pipeline is designed to analyze a FASTQ file, build a De Bruijn graph, and use a Bloom filter to track k-mers efficiently. Here's a breakdown of how it works:
### 1. Rust Program (main.rs)

This program processes a FASTQ file, extracts k-mers, constructs a De Bruijn graph, and builds a Bloom filter.
Key Components

    De Bruijn Graph Representation
        Uses k-mers (subsequences of length k) as nodes.
        Each node stores edges to next k-mers found in the read sequences.

    Bloom Filter Implementation
        A probabilistic data structure that allows fast membership checks with minimal memory usage.
        Uses multiple hash functions to reduce false positives.
        Supports parallel insertion of k-mers.

### Execution Flow

    Parse CLI arguments (FASTQ path, k-mer size, output directory).
    Read the FASTQ file using needletail (a Rust bioinformatics library).
    Extract k-mers from the sequences and construct the De Bruijn graph.
    Insert k-mers into the Bloom filter.
    Check if a specific k-mer (ACGT) is in the Bloom filter.
    Save results as JSON files:
        graph.json: Contains the De Bruijn graph.
        bloom.json: Contains the Bloom filter.

### 2. Nextflow Pipeline (main.nf)

This automates the execution of the Rust program in a workflow.
Processes in the Pipeline

    Compile the Rust Code (compile process)
        Uses cargo build --release to generate an optimized Rust binary.
        Saves the binary (debruijn_bloom).

    Run the Analysis (analysis process)
        Runs the compiled Rust binary with the provided FASTQ file.
        Produces graph.json and bloom.json.

### Workflow Execution

    The workflow section first compiles the Rust code, then runs analysis on the compiled binary.

### 3. Output Explanation (output.txt)

Reading from FASTQ: example.fastq, k-mer=31
Loaded 1 reads. Building de Bruijn graph & Bloom filter...
Bloom filter contains 'ACGT'? false
De Bruijn graph saved to results/graph.json
Bloom filter saved to results/bloom.json.

    Reads example.fastq with k=31.
    Constructs the De Bruijn graph and Bloom filter.
    Checks if "ACGT" is in the Bloom filter (not found).
    Saves results in results/graph.json and results/bloom.json.

### Key Takeaways

    De Bruijn graphs are useful for genome assembly and sequence analysis.
    Bloom filters efficiently store k-mers with low memory usage.
    Parallel processing (rayon) speeds up k-mer extraction and Bloom filter insertion.
    Nextflow automates the pipeline, making it reproducible and scalable.
