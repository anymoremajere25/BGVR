### 5.4. De Novo Assembly Approaches  
**experiment_5_4_1**

This Rust program illustrates a chunk-based method for constructing a k-mer count table and generating a minimal de Bruijn graph from large FASTQ inputs. By dividing the FASTQ reading process into manageable chunks, it avoids overloading memory. Each chunk is processed concurrently using Rayon, and the partial de Bruijn graphs are serialized to disk. In high-performance computing (HPC) or cloud-based environments, this chunking step can run in parallel across different portions of the input, producing separate partial graphs that can later be merged.

Upon parsing command-line arguments with `clap`, the program creates an output directory to store the partial de Bruijn graphs. It employs the `parse_fastx_file` function from `needletail` to stream through the FASTQ records, reading a fixed number of records (defined by `chunk_size`) at a time. This chunk is then processed in parallel using Rayon’s `.par_iter()`, distributing the task of k-mer counting across available CPU cores.

During the processing of each chunk, the program constructs a local `FnvHashMap` to store k-mer counts and uses it to build a partial de Bruijn graph, treating each k-mer’s prefix as a node and the final base as an edge. This partial graph is serialized to disk using `bincode`. After all chunks are processed, the program reads and merges the partial graphs into one complete de Bruijn graph. A final thresholding pass removes edges with counts below a specified coverage value, resulting in a cleaner graph. This approach of generating partial outputs followed by a merging step is commonly used in HPC pipelines, enabling the efficient scaling to large datasets while maintaining safe parallel processing.

**File Contents:**

    experiment_5_4_1/
        Cargo.toml (dependency file)  
        experiment_5_4_1/src/
            main.rs (Rust script)  
            reads.fq (FASTQ file)  
            reference.fa (FASTA file)


**Dependencies:**

```toml
[dependencies]
anyhow = "1.0"
rayon = "1.8"
needletail = "0.6.3"
fnv = "1.0"
serde = { version = "1.0", features = ["derive"] }
serde_json = "1.0"
bincode = "2.0.1"
clap = { version = "4.4", features = ["derive"] }
```

### Explanation of Output  

**Step 1: Reading the FASTQ File**  
The program reads the FASTQ file in chunks (default: 10,000 sequences per chunk). Each sequence is processed to extract k-mers (default: k=31).

**Step 2: Counting k-mers**  
Each chunk is processed in parallel using Rayon, where:
- Every sequence is scanned for overlapping k-mers.
- A hashmap stores the count of each k-mer.
- Only k-mers that meet a minimum threshold (default: 2 occurrences) are included in the graph.

**Step 3: Constructing Partial de Bruijn Graphs**  
A partial de Bruijn graph is constructed for each chunk and saved to disk. These graphs are stored in binary format using `bincode` inside the directory `partial_kmer_maps/`.

**Step 4: Merging Partial Graphs**  
The program reads the partial graphs and merges them into a final graph. Any edges (k-1-mers → next base) with counts below the threshold are removed.

**Step 5: Writing Final de Bruijn Graph**  
The final de Bruijn graph is stored as a binary file (`final_debruijn.bin`).

### Example Output

Assuming the following command is run:

```bash
cargo run -- --input example.fastq --k 5 --threshold 2 --chunk_size 5000
```

**Output:**

```
Processed chunk 0 with 5000 records, wrote partial de Bruijn to "partial_kmer_maps/partial_debruijn_0.bin"
Processed chunk 1 with 5000 records, wrote partial de Bruijn to "partial_kmer_maps/partial_debruijn_1.bin"
Processed chunk 2 with 4200 records, wrote partial de Bruijn to "partial_kmer_maps/partial_debruijn_2.bin"
Merging partial de Bruijn graphs...
Final de Bruijn graph has 1,234 prefix nodes. Written to "final_debruijn.bin".
```

**Generated Files:**

- `partial_kmer_maps/partial_debruijn_0.bin`
- `partial_kmer_maps/partial_debruijn_1.bin`
- `partial_kmer_maps/partial_debruijn_2.bin`
- `final_debruijn.bin`

### Interpretation of Output

- **Processed chunk X with Y records:** Each chunk of the FASTQ file is processed and a partial graph is saved.
- **Final de Bruijn graph has N prefix nodes:** The final graph contains N unique (k-1)-mer prefixes, which indicates N unique k-mers were found.
- **final_debruijn.bin:** The merged and filtered de Bruijn graph, saved in binary format for future analysis.

### Conclusion

The program efficiently constructs a de Bruijn graph in parallel and with minimal memory usage. The chunking mechanism allows large FASTQ files to be processed without high memory demands. The final graph can be used for genome assembly, sequence error correction, or other bioinformatics analyses.
