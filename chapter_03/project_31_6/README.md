
## 3.1. Introduction to Data Processing in Bioinformatics
### project_31_6

This Rust program processes a FASTA file by reading each sequence record, computing its GC content, and saving the results in JSON format. The program efficiently parses sequences using the needletail crate, which allows for fast and memory-efficient processing. Each sequence is converted into a structured format containing:

    Sequence ID
    Sequence Length
    GC Content (percentage of G and C nucleotides)

Once processed, the data is written to an output JSON file for further analysis.
### Algorithm and Complexity Analysis

The GC content for each sequence is computed by counting occurrences of G and C, then dividing by the sequence length. Given:

    nn = number of sequences
    ℓℓ = average length of a sequence

The GC content calculation runs in O(\ell) for each sequence, leading to an overall complexity of O(n × \ell). This makes the approach highly efficient for large genomic datasets.

If parallelization were introduced, the GC computation could be distributed across multiple CPU cores, reducing execution time without changing the asymptotic complexity.
Project Structure

### project_31_6/
    Cargo.toml         # Rust package configuration and dependencies  
    src/
        main.rs       # Main Rust script  
        example.fasta # Sample FASTA file containing DNA sequences  
        output.txt    # Console output  
        output.json   # Final processed results in JSON format  

### How to Run the Program

In Cursor , use the following command to execute and save the output:

cargo run > output.txt

This runs main.rs and stores the output in output.txt.
Required Dependencies (Cargo.toml)

[dependencies]
needletail = "0.6.3"  # Efficient FASTA/FASTQ parser  
serde = { version = "1.0", features = ["derive"] }  # For struct serialization  
serde_json = "1.0"  # For writing JSON output  

### Detailed Breakdown of Program Execution
1. Reading the FASTA File

    The program expects a FASTA file as input.
    If no file path is provided, it defaults to "example.fasta".
    Uses needletail::parse_fastx_file to efficiently read and process sequences.

2. Extracting Sequence Information

For each sequence, the program computes:

    GC content: The fraction of G and C bases in the sequence.
    Sequence length: Total number of nucleotides.
    Sequence ID: Extracted from the FASTA header.

This data is stored in the SeqRecord struct:

struct SeqRecord {
    id: String,
    seq_len: usize,
    gc_content: f64,
}

3. GC Content Calculation

The function calc_gc_content computes GC content by filtering and counting G and C nucleotides:

fn calc_gc_content(seq: &[u8]) -> f64 {
    let gc_count = seq.iter().filter(|&&c| c == b'G' || c == b'C').count();
    gc_count as f64 / (seq.len().max(1) as f64)
}

4. Writing Results to JSON

    Once all sequences are processed, the results are serialized into output.json using serde_json.
    The program also prints the number of processed sequences to the console.

### Example Input and Output
Input: example.fasta

>seq1
AGCTGCGCGGATCAGCGT
>seq2
CGTAGCTAGGCTAGGCTA

Console Output: output.txt

Processed 2 sequences

JSON Output: output.json

[
  {"id": "seq1", "seq_len": 18, "gc_content": 0.6666666666666666},
  {"id": "seq2", "seq_len": 18, "gc_content": 0.5555555555555556}
]

### Conclusion

    The program successfully reads, processes, and analyzes sequences from a FASTA file.
    It computes sequence length and GC content and saves the results in JSON format.
    The approach is efficient (O(n × \ell)) and scalable for large datasets.
    The output confirms successful parsing, with two sequences processed.

Key Takeaways

✔ Efficient FASTA parsing using needletail
✔ GC content analysis for genomic sequences
✔ JSON output for structured data storage
✔ Scalable approach for large bioinformatics datasets

This implementation can be extended for parallel processing, large-scale genomic studies, or integration with bioinformatics pipelines.





