**3.2. Sequence Data Structures and String Algorithms**  
**Experiment_3_2_3**  

### 1. Nextflow
The Nextflow script below defines two processes: one for compiling a Rust binary (optional if using a precompiled artifact) and another for executing a partial suffix array program on a given FASTA file. Each segment of the FASTA sequence is processed in parallel using Rayon, generating multiple partial arrays, which are serialized as JSON. In large-scale applications, ephemeral nodes may handle individual sequence chunks before merging the outputs into a unified suffix array. This Rust-based approach demonstrates chunking, parallel sorting, and offset adjustments efficiently.

### 2. Rust Implementation
The Rust code processes the FASTA file by first removing header lines (lines starting with '>') and concatenating the remaining sequence into a single string. The sequence is then split into manageable chunks (e.g., 1 million bases), utilizing `as_bytes().chunks(chunk_size)`. For each chunk, a naive suffix array is created by enumerating all possible substring start positions and sorting them. While this method is impractical for extremely large datasets, it demonstrates the core concept effectively. The parallel processing is facilitated by Rayon’s `.par_iter()`, distributing workload across CPU cores for efficiency. The output is serialized into `partial_suffix_arrays.json`, which can be merged or analyzed further.

To execute the code locally, place `main.nf` (Nextflow script) in the same directory as the Rust project, then run:
```powershell
nextflow run main.nf --fasta /path/to/large_sequence.fa
```
This triggers the compilation of the Rust binary (if necessary) and runs the analysis on the provided FASTA file. The output, `partial_suffix_arrays.json`, contains the serialized suffix arrays.

For large-scale environments, each job may produce its own JSON file, later merged by Nextflow or an external process. Boundary-spanning operations like full BWT construction or k-mer overlap handling can be incorporated by including extra bases at chunk boundaries. Given that the naive approach sorts all suffix positions in memory, it is resource-intensive. Production-ready solutions often use advanced algorithms like SA-IS, chunk partitioning, or distributed indexing frameworks. Robust error handling, file validation, and resource management are essential for processing large datasets efficiently. Nextflow orchestrates this workflow across local or cloud-based HPC environments, while Rust ensures thread-safe parallel processing.

### Project Directory Structure
```
experiment_3_2_3/
    Cargo.toml  (Dependencies configuration)
    src/
        main.rs  (Rust script)
        main.nf  (Nextflow script)
    large_sequence.fa  (FASTA input file)
    partial_suffix_arrays.json  (JSON output file)
    output.txt  (Execution log)
```

### Execution Instructions
Run the following command in PowerShell:
```powershell
cargo run main.nf | tee output.txt
```
This executes `main.nf` and saves the output to `output.txt`.

#### Dependencies
```toml
[dependencies]
rayon = "1.7"
serde = { version = "1", features = ["derive"] }
serde_json = "1"
```

### Output Explanation
#### 1. `output.txt`
This log file contains messages indicating the number of partial suffix arrays generated and stored in `partial_suffix_arrays.json`. Example output:
```
Generated 1 partial array in partial_suffix_arrays.json
```
This suggests the FASTA file was processed in a single chunk, likely due to its size fitting within the predefined 1,000,000-character limit.

#### 2. `partial_suffix_arrays.json`
This JSON file stores serialized partial suffix arrays computed by the Rust program. Example snippet:
```json
[
  {
    "start_pos": 0,
    "suffix_positions": [589, 979, 692, 1126, 1045, 365, 718, 54, 104, ...]
  }
]
```
- **start_pos**: Denotes the chunk's starting position in the sequence.
- **suffix_positions**: A sorted list of suffix start positions for this chunk. For example, `589` indicates that the suffix starting at index `589` in the sequence is the lexicographically smallest in this chunk.

### Execution Breakdown
1. **FASTA Processing**: The program reads the file, removes headers, and extracts the sequence.
2. **Chunking**: The sequence is divided into 1,000,000-character segments.
3. **Parallel Suffix Array Generation**:
   - Enumerates suffix start positions for each chunk.
   - Sorts suffixes lexicographically.
   - Uses Rayon for parallel processing across CPU cores.
4. **Nextflow Automation**:
   - Automates compilation and execution of Rust code.
   - Saves the resulting suffix arrays in JSON format for further processing.

### Key Insights
1. **Successful Execution**:
   - The program processed the FASTA file and generated 1 partial suffix array, as confirmed by `output.txt` and `partial_suffix_arrays.json`.
2. **Parallel Efficiency**:
   - Rayon efficiently distributes workload across CPU cores, which scales well for large datasets.
3. **Suffix Array Construction**:
   - Outputs sorted indices of suffix positions, forming the basis for full suffix array merging.
4. **Nextflow Workflow Management**:
   - Seamlessly handles the compilation and execution of Rust-based computations.
5. **Scalability**:
   - Designed to handle large FASTA files efficiently. With increasing sequence sizes, parallelization becomes even more beneficial.

### Conclusion
This pipeline successfully demonstrates a parallelized suffix array generation method using Rust and Nextflow. The approach is scalable and provides a foundation for processing much larger biological datasets efficiently.

