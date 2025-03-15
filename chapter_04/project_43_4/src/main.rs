use rayon::prelude::*;

/// A configurable pattern describing the TATA box or a variant. Supports mismatch tolerance.
#[derive(Debug)]
pub struct TATAPattern {
    /// The core consensus we are looking for. E.g., "TATA(A/T)A" → "TATA?A" (where '?' can be A/T).
    consensus: Vec<Vec<char>>,
    /// Maximum number of mismatches allowed, e.g., 0 means exact matching only.
    pub max_mismatches: usize,
}

impl TATAPattern {
    /// Construct a default TATA pattern with zero mismatches and the exact consensus "TATA(A|T)A".
    pub fn default_tata() -> Self {
        // For a 6-base motif:
        //   T, A, T, A, (A|T), A
        // We store each position as a vector of acceptable chars.
        let consensus = vec![
            vec!['T'],       // position 0
            vec!['A'],       // position 1
            vec!['T'],       // position 2
            vec!['A'],       // position 3
            vec!['A', 'T'],  // position 4 can be A or T
            vec!['A'],       // position 5
        ];
        TATAPattern {
            consensus,
            max_mismatches: 0,
        }
    }

    /// Construct a TATAPattern from an arbitrary set of acceptable nucleotides at each position.
    /// For example, if you want more flexibility in the 3rd position, supply multiple chars for that position.
    pub fn new(consensus: Vec<Vec<char>>, max_mismatches: usize) -> Self {
        TATAPattern {
            consensus,
            max_mismatches,
        }
    }

    /// Returns the length of the pattern in nucleotides.
    pub fn len(&self) -> usize {
        self.consensus.len()
    }

    /// Checks whether `window` matches this pattern within `max_mismatches` tolerance.
    /// Case-insensitive: 'T' == 't', etc.
    fn matches_window(&self, window: &[char]) -> bool {
        let mut mismatches = 0;
        for (i, &c) in window.iter().enumerate() {
            // Convert to uppercase for case-insensitive comparison
            let c_up = c.to_ascii_uppercase();
            if !self.consensus[i].contains(&c_up) {
                mismatches += 1;
                if mismatches > self.max_mismatches {
                    return false;
                }
            }
        }
        true
    }
}

/// Find all matches of a TATA-like pattern within a single sequence.
pub fn find_tata_boxes(sequence: &str, pattern: &TATAPattern) -> Vec<usize> {
    let mut positions = Vec::new();
    let needed_len = pattern.len();
    if needed_len == 0 {
        // If the pattern length is zero, return nothing (or handle as error).
        return positions;
    }
    let chars: Vec<char> = sequence.chars().collect();
    if chars.len() < needed_len {
        return positions;
    }
    for start in 0..=(chars.len() - needed_len) {
        let window = &chars[start..start + needed_len];
        if pattern.matches_window(window) {
            positions.push(start);
        }
    }
    positions
}

/// Scan multiple DNA sequences in parallel, returning a vector of match positions for each sequence.
pub fn find_tata_boxes_parallel(sequences: &[String], pattern: &TATAPattern) -> Vec<Vec<usize>> {
    // Use rayon parallel iterator to process each sequence concurrently.
    sequences
        .par_iter()
        .map(|seq| find_tata_boxes(seq, pattern))
        .collect()
}

/// Demonstration main function, showing usage on multiple sequences.
fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Example sequences
    let seqs = vec![
        "GGTTTATATAAACTATAATTTTACGT".to_string(),
        "tttatacccggttttataAa".to_string(),
        "AAAAATATA".to_string(),
        "NoTATAhere".to_string(),
    ];

    // Example 1: default TATA pattern, no mismatches
    let default_pattern = TATAPattern::default_tata();

    // Example 2: custom pattern with mismatch tolerance
    // For demonstration, let's allow a mismatch in the 4th position.
    // Pattern: T A T A (A/T), A, but 1 mismatch allowed across these 6 positions.
    let custom_consensus = vec![
        vec!['T'],
        vec!['A'],
        vec!['T'],
        vec!['A'],
        vec!['A', 'T'],
        vec!['A'],
    ];
    let custom_pattern = TATAPattern::new(custom_consensus, 1);

    // Run the parallel search for default pattern
    let all_matches_default = find_tata_boxes_parallel(&seqs, &default_pattern);
    // Run the parallel search for custom pattern
    let all_matches_custom = find_tata_boxes_parallel(&seqs, &custom_pattern);

    println!("Sequences: {:#?}\n", seqs);

    for (i, positions) in all_matches_default.iter().enumerate() {
        println!("Seq {} - TATA default pattern matches: {:?}", i, positions);
    }
    println!();

    for (i, positions) in all_matches_custom.iter().enumerate() {
        println!("Seq {} - TATA custom pattern (1 mismatch) matches: {:?}", i, positions);
    }

    Ok(())
}
