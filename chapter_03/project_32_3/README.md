Nextflow Pipeline (main.nf)

This script defines an automated pipeline in Nextflow, handling:

    Compiling a Rust program (compile process)
    Running the compiled Rust binary on a FASTA file (analysis process)

Pipeline Overview

    Inputs
        A FASTA file (large_sequence.fa) containing a DNA sequence.
    Outputs
        A JSON file (partial_suffix_arrays.json) storing a partial suffix array.

🚀 Process Breakdown
🔹 compile Process

    Compiles the Rust program using cargo build --release
    Copies the compiled binary (project_32_3) to partial_sa for easier execution.
    Ensures it runs once using Channel.of(1).

🔹 analysis Process

    Runs the compiled binary (./partial_sa) on the provided FASTA file.
    Generates a partial suffix array and saves it as partial_suffix_arrays.json.

🔹 workflow Execution Order

    compile() → Builds the Rust binary.
    analysis() → Runs the binary on the FASTA file.

2️⃣ Rust Code (main.rs)

This Rust program reads a DNA sequence, processes it in parallel, and generates a partial suffix array.
🚀 Step-by-Step Execution
🔹 Step 1: Read the FASTA File

    Reads the input FASTA file (large_sequence.fa).
    Ignores header lines (lines starting with >).
    Extracts only the DNA sequence.

🔹 Step 2: Chunk the Sequence

    Splits the DNA sequence into 1 million base pair chunks for parallel processing.

🔹 Step 3: Compute Partial Suffix Arrays in Parallel

    Uses Rayon for parallel processing.
    Each chunk generates a sorted list of suffix positions.
    The starting position of each chunk is tracked.

🔹 Step 4: Serialize and Save Output

    Converts the partial suffix array into JSON format.
    Saves it as partial_suffix_arrays.json.

3️⃣ Expected Output
Example Input (FASTA File)

>example_sequence
ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGC

Example Output (partial_suffix_arrays.json)

[
  {
    "start_pos": 0,
    "suffix_positions": [5, 10, 15, 0, 20, 25]  
  }
]

    start_pos: The chunk's starting position in the original sequence.
    suffix_positions: The sorted suffix positions within the chunk.

4️⃣ Explanation of Generated Message

Generated 1 partial arrays in partial_suffix_arrays.json

    1 partial suffix array was created from 1 chunk of DNA.
    The data is saved in JSON format for further analysis.

5️⃣ Summary

    ✅ Nextflow automates compilation and execution.
    ✅ Rust program computes partial suffix arrays in parallel.
    ✅ JSON output stores suffix positions for later merging.
