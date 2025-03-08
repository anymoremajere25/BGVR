use rayon::prelude::*;
use serde::{Serialize, Deserialize};
use std::fs;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::path::PathBuf;
use std::error::Error;

/// Represents a partial suffix array for a chunk of the sequence.
#[derive(Debug, Serialize, Deserialize)]
struct PartialSuffixArray {
    start_pos: usize,
    suffix_positions: Vec<usize>,
}

fn main() -> Result<(), Box<dyn Error>> {
    // 1) Take FASTA path from CLI args (or default to "large_sequence.fa").
    let fasta_path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "large_sequence.fa".to_string());
    // 2) Read contents, skipping header lines in typical FASTA format.
    let contents = fs::read_to_string(PathBuf::from(&fasta_path))?;
    let sequence: String = contents
        .lines()
        .filter(|line| !line.starts_with('>'))
        .collect();

    // 3) Define chunk size and partition the sequence into manageable slices.
    let chunk_size = 1_000_000;
    let chunks: Vec<_> = sequence
        .as_bytes()
        .chunks(chunk_size)
        .map(|c| String::from_utf8_lossy(c).to_string())
        .collect();

    // For demonstration, we track how many partial arrays we build in total.
    let partial_count = AtomicUsize::new(0);

    // 4) Build a partial suffix array for each chunk in parallel.
    let partial_arrays: Vec<PartialSuffixArray> = chunks
        .par_iter()
        .enumerate()
        .map(|(idx, chunk)| {
            let offset = idx * chunk_size;
            // Create an array of suffix starting positions, then sort by substring.
            let mut suffix_positions: Vec<usize> = (0..chunk.len()).collect();
            suffix_positions.sort_by_key(|&pos| &chunk[pos..]);
            partial_count.fetch_add(1, Ordering::Relaxed);

            PartialSuffixArray {
                start_pos: offset,
                suffix_positions,
            }
        })
        .collect();

    // 5) Serialize partial arrays to JSON for potential merging later.
    let serialized = serde_json::to_string_pretty(&partial_arrays)?;
    fs::write("partial_suffix_arrays.json", serialized)?;

    println!(
        "Generated {} partial arrays in partial_suffix_arrays.json",
        partial_count.load(Ordering::Relaxed)
    );
    Ok(())
}
