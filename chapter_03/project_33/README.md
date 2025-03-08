# 3.3. Graph Data Structures for Genome Assembly and Beyond

### project_33

This project demonstrates how to build partial de Bruijn graphs from FASTQ data in Rust. The implementation utilizes needletail for efficient sequence parsing and Petgraph to represent overlapping k-mers in an undirected graph. By dividing reads into chunks and processing them in parallel using Rayon, each segment contributes to a partial de Bruijn graph, which is serialized for further merging or analysis. This design efficiently handles large genomic datasets by distributing computational workload across multiple cores while optimizing memory usage.

### Overview

#### The program:

Reads all sequences from a FASTQ file.

Divides them into equal-sized chunks.

Processes each chunk in parallel by extracting k-mers from reads.

Constructs a graph, adding nodes for k-mers and edges between consecutive k-mers.

Serializes and saves partial de Bruijn graphs to a JSON file.

This modular approach enables intermediate results to be stored and later merged into a complete de Bruijn graph if needed.

### Project Structure

project_33/
│── Cargo.toml                   # Dependencies configuration
│── src/
│   ├── main.rs                  # Rust script for graph construction
│   ├── reads.fq.rar              # Compressed input FASTQ file
│   ├── partial_debruijn_graphs.json.rar  # Compressed output JSON file
│   ├── output.txt                # Execution log file

### How to Run

Run the following command in cursor terminal:

cargo run > output.txt

This executes main.rs and saves the program's output in output.txt.

### Dependencies (Cargo.toml)

[dependencies]
rayon = "1.7"
needletail = "0.6"
serde = { version = "1", features = ["derive"] }
serde_json = "1"
petgraph = "0.7.1"

### Output Explanation

1. output.txt

Wrote partial de Bruijn graphs to partial_debruijn_graphs.json

This message confirms that the graph construction completed successfully and was saved to partial_debruijn_graphs.json.

2. partial_debruijn_graphs.json

This JSON file contains the constructed partial de Bruijn graphs with the following structure:

nodes: List of unique k-mers (e.g., 31-mers).

edges: List of tuples (i, j) representing edges between k-mers.

k: Length of each k-mer (e.g., 31).

### Example JSON Output:

{
  "nodes": [
    "ACGGAGGATGCGAGCGTTATCCGGATTTATT",
    "CGGAGGATGCGAGCGTTATCCGGATTTATTG",
    "GGAGGATGCGAGCGTTATCCGGATTTATTGG",
    "GAGGATGCGAGCGTTATCCGGATTTATTGGG",
    "AGGATGCGAGCGTTATCCGGATTTATTGGGT",
    "GGATGCGAGCGTTATCCGGATTTATTGGGTT"
  ],
  "edges": [
    [0, 1],
    [1, 2],
    [2, 3],
    [3, 4],
    [4, 5]
  ],
  "k": 31
}

3. Algorithm Breakdown

Step 1: Read Input FASTQ/FASTA File

The program reads sequencing reads from reads.fq.

It uses the needletail library for efficient FASTA/FASTQ parsing.

Step 2: Partition Reads into Chunks

Reads are split into smaller chunks (100,000 reads per chunk).

This enables efficient parallel processing.

Step 3: Construct Partial de Bruijn Graphs

Each chunk undergoes:

k-mer extraction (substring of length k = 31).

Node addition to the graph.

Edge creation between consecutive k-mers.

Resulting graph storage.

Step 4: Parallel Processing

Rayon is used to process multiple chunks in parallel.

Each chunk contributes to a partial de Bruijn graph.

Step 5: Output Partial Graphs

Partial graphs are serialized to JSON (partial_debruijn_graphs.json).

This enables further processing, such as merging into a global graph.

### Conclusion

1. Successful Construction of Partial de Bruijn Graphs

The output JSON confirms successful graph generation.

2. Efficient Parallelization

Using Rayon allows for parallel processing, improving performance for large datasets.

3. Scalability

The chunk-based approach efficiently handles massive sequencing datasets without exhausting memory.

4. Accurate Graph Representation

Nodes represent k-mers, and edges properly reflect overlapping k-mers.

5. Next Steps

Merge partial graphs to form a complete de Bruijn graph.

Apply results to genome assembly, error correction, or variant detection.

### Final Thoughts

This pipeline efficiently constructs partial de Bruijn graphs from sequencing reads.

Parallel execution improves performance and scales well for large sequencing datasets.

The graph-based approach is useful for genome assembly and sequence analysis.


