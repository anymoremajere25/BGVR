## 5.3. Alignment and Mapping Algorithms
### Experiment 5.3.1 {fmindex aligner}

This Rust program implements an FM-index-based read alignment tool. Here's how it works:

### **1. Argument Parsing (`Args` struct)**
- The program uses `clap` to parse command-line arguments for:
  - `reference`: The reference genome file.
  - `reads`: The file containing reads to align.
  - `index_outdir`: Directory to store FM-index chunks.
  - `reference_chunk_size`: Chunk size for reference genome (default: 50).
  - `fm_sampling_rate`: FM-index sampling rate (default: 3).
  - `alignment_output`: Output file for alignments.

### **2. Preprocessing the Reference Genome**
- Converts reference genome to uppercase.
- Removes invalid characters (keeps `A, C, G, T, N`).
- Ensures it ends with a `$` sentinel character.

### **3. FM-Index Construction**
- **Suffix Array (`suffix_array`)**: Constructs the suffix array of the reference genome.
- **Burrows-Wheeler Transform (`bwt`)**: Computes the BWT using the suffix array.
- **Less & Occurrence (`Less`, `Occ`)**: Constructs auxiliary data structures needed for FM-index.
- **FM-Index (`FMIndex`)**: Uses `bio` crate to construct the FM-index.

### **4. Read Alignment Using FM-Index**
- Reads sequences are processed in parallel using `rayon`.
- Each read is searched against the FM-index using **backward search**.
- If a match is found, the corresponding positions in the reference genome (from suffix array) are recorded.
- If no match is found, an empty result is stored.

### **5. Output Alignment Results**
- The results are serialized into a JSON file using `serde_json`.

### **6. Dependencies (`Cargo.toml`)**
- `anyhow`: Error handling.
- `rayon`: Parallel computation.
- `bio`: Provides FM-index, BWT, and suffix array functionality.
- `serde` & `serde_json`: Serialization.
- `clap`: Command-line argument parsing.

### **Expected Outcome**
- The program reads a reference genome and a set of reads.
- It builds an FM-index and aligns reads efficiently.
- The alignment results are stored in a JSON file.

Would you like a test example or further clarifications? 😊
