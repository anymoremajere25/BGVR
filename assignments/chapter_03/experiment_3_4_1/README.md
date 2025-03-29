### 3.4. Indexing and Searching in Large-Scale Biological Datasets
**experiment_3_4_1**

This implementation demonstrates a minimal Bloom filter in Rust for k-mers extracted from FASTA or FASTQ files. The program reads an input file, splits the sequences into k-mers in parallel using Rayon, and sequentially inserts them into a Bloom filter. The filter is then serialized to JSON for further analysis. The parallel k-mer extraction efficiently distributes computation across CPU cores, while single-threaded Bloom filter insertion ensures thread safety.

After processing all sequences, overlapping k-mers are generated in a thread-safe manner and inserted into a Bloom filter. The filter, configured based on a target false-positive rate, utilizes a bit vector and a multi-seed hashing strategy for efficient storage and retrieval. The serialized Bloom filter allows for querying k-mer membership and can be merged with partial filters from other sources.

For high-performance computing (HPC) environments, ephemeral jobs can generate partial Bloom filters that are later combined using bitwise operations. This facilitates fast membership queries in large-scale genomic applications, such as read classification or reference indexing. Careful tuning of filter parameters (number of bits and hash functions) balances memory efficiency and false-positive rates. In large-scale applications, prefix-based k-mer distribution across nodes can optimize load balancing.

### Project Structure
```
experiment_3_4_1/
    Cargo.toml  (Dependencies and build configuration)
    src/
        main.rs  (Rust source code)
    reads.fq.rar  (Compressed sequencing reads)
    bloom.json  (Serialized Bloom filter)
    output.txt  (Execution output)
```

### Running the Code
Run the following command in PowerShell:
```sh
cargo run | tee output.txt
```
(This executes `main.rs` and saves the output to `output.txt`.)

### Dependencies
```toml
[dependencies]
rayon = "1.7"
needletail = "0.6"
serde = { version = "1", features = ["derive"] }
serde_json = "1"
bitvec = "1.0.1"
```

### Explanation of the Output

#### `output.txt`
```
Bloom contains first k-mer 'ACGGAGGATGCGAGCGTTATCCGGATTTATT': true
Constructed Bloom filter (10,650,458 k-mers processed), result written to bloom.json
```

**Breakdown:**
- **10,650,458 k-mers** were processed and stored in the Bloom filter.
- The first extracted k-mer (`ACGGAGGATGCGAGCGTTATCCGGATTTATT`) is confirmed to be present.
- The Bloom filter is serialized and saved in `bloom.json`.

#### `bloom.json`
Stores the Bloom filter configuration and data:
```json
{
  "bits": [4, 0, 48, 128, 10, 0, 0, 66, ...],
  "num_bits": 10000000,
  "num_hashes": 3,
  "k": 31
}
```
- **`bits`**: The bit array representing stored k-mers.
- **`num_bits`**: Total bit space allocated (10 million bits).
- **`num_hashes`**: Number of hash functions applied per k-mer (3).
- **`k`**: Length of each k-mer (31 nucleotides).

### Process Breakdown
1. **Read Sequencing Data**
   - Reads sequences from `reads.fq`.
   - Extracts all overlapping 31-mers.

2. **Initialize Bloom Filter**
   - Allocates a 10MB bit vector.
   - Uses three hash functions for insertion.

3. **Parallel k-mer Extraction**
   - Uses `rayon` to speed up processing.
   - Extracts 10,650,458 k-mers.

4. **Insert k-mers into Bloom Filter**
   - Each k-mer is hashed and mapped to multiple bit positions.

5. **Check k-mer Membership**
   - Confirms existence of an example k-mer.

6. **Serialize and Save Bloom Filter**
   - Outputs filter to `bloom.json`.

### Key Insights
- **Memory Efficiency:** Stores 10M+ k-mers in just 10MB.
- **Probabilistic Nature:** Allows false positives but guarantees no false negatives.
- **Parallel Processing:** Leverages multiple CPU cores for faster extraction.
- **Scalability:** Can be extended for genome-scale datasets.

### Next Steps
- Use the Bloom filter for **genome assembly** and **error filtering**.
- Optimize hash functions and bit space for larger datasets.
- Implement **distributed Bloom filters** for even greater scalability.

### Final Thoughts
This implementation provides a compact and efficient k-mer membership test using a Bloom filter. The combination of Rust, parallel processing, and Nextflow workflows enables handling large genomic datasets with high performance.


