## 4.3 Motif Discovery and regulatory element identification

### Explanation [project_43_4] Rust Program

This Rust program is designed to find occurrences of a TATA-like motif in DNA sequences, both with exact matches and allowing one mismatch. It uses parallel processing (via Rayon) to efficiently search across multiple sequences.
#### 1. Overview of How It Works
Step 1: Define the TATA Pattern (TATAPattern struct)

    The struct stores a consensus sequence (the expected motif).

    It also allows a configurable number of mismatches.

Key Methods in TATAPattern

    default_tata(): Defines a strict TATA(A/T)A consensus with zero mismatches.

    new(): Allows defining a flexible TATA variant with custom mismatches.

    matches_window(): Checks if a substring (window) matches the consensus sequence, counting mismatches.

Step 2: Find Motifs in Sequences

    find_tata_boxes(): Scans a single DNA sequence for motif matches.

    find_tata_boxes_parallel(): Uses Rayon to scan multiple sequences concurrently.

Step 3: Run the Program (main.rs)

    Load Example DNA Sequences

let seqs = vec![
    "GGTTTATATAAACTATAATTTTACGT".to_string(),
    "tttatacccggttttataAa".to_string(),
    "AAAAATATA".to_string(),
    "NoTATAhere".to_string(),
];

Define the TATA motif search patterns

    Exact match pattern (default).

    Custom pattern allowing 1 mismatch.

Search for matches in parallel

let all_matches_default = find_tata_boxes_parallel(&seqs, &default_pattern);
let all_matches_custom = find_tata_boxes_parallel(&seqs, &custom_pattern);

Print results

    println!("Seq {} - TATA default pattern matches: {:?}", i, positions);
    println!("Seq {} - TATA custom pattern (1 mismatch) matches: {:?}", i, positions);

### 2. Output Analysis

The program prints:
Seq 0 - TATA default pattern matches: [4, 6]
Seq 1 - TATA default pattern matches: [14]
Seq 2 - TATA default pattern matches: []
Seq 3 - TATA default pattern matches: []

Seq 0 - TATA custom pattern (1 mismatch) matches: [2, 4, 6, 13]
Seq 1 - TATA custom pattern (1 mismatch) matches: [0, 12, 14]
Seq 2 - TATA custom pattern (1 mismatch) matches: [3]
Seq 3 - TATA custom pattern (1 mismatch) matches: []


Seq 0 - TATA default pattern matches: [4, 6]
Seq 1 - TATA default pattern matches: [14]
Seq 2 - TATA default pattern matches: []
Seq 3 - TATA default pattern matches: []

Seq 0 - TATA custom pattern (1 mismatch) matches: [2, 4, 6, 13]
Seq 1 - TATA custom pattern (1 mismatch) matches: [0, 12, 14]
Seq 2 - TATA custom pattern (1 mismatch) matches: [3]
Seq 3 - TATA custom pattern (1 mismatch) matches: []

This means:

    The exact "TATA" pattern is found in Seq 0 (positions 4,6) and Seq 1 (position 14).

    The 1-mismatch version finds more matches across the sequences.

### 3. Optimization & Next Steps
✅ Strengths

    Efficient: Uses parallel computing (Rayon).

    Flexible: Can customize mismatches and pattern variability.

    Case-insensitive motif matching.

⚡ Possible Improvements

    Allow configurable mismatch tolerance (not just 0 or 1).

    Use a Position Weight Matrix (PWM) instead of strict character matching.

    Enhance sequence input handling (e.g., load from a FASTA file).
