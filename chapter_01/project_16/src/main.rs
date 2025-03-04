use rayon::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader, Result};

/// A helper function to load a large genomic sequence from a file into a String.
/// In real use, you might stream data or handle multi-GB files in chunks.
fn load_genome(path: &str) -> Result<String> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    // Collect all lines, skipping any FASTA headers (lines starting with '>').
    let mut genome = String::new();
    for line in reader.lines() {
        let line = line?;
        if !line.starts_with('>') {
            genome.push_str(&line.trim());
        }
    }
    Ok(genome)
}

/// Constructs a naive suffix array for a given &str by sorting suffixes.
/// Returns a vector of starting indices in sorted order.
fn build_suffix_array(seq: &str) -> Vec<usize> {
    let mut suffixes: Vec<usize> = (0..seq.len()).collect();
    // Sort by comparing slices directly; O(n log n) with expensive comparisons.
    suffixes.sort_by(|&a, &b| seq[a..].cmp(&seq[b..]));
    suffixes
}

/// Splits the genome into `num_chunks` pieces, builds partial suffix arrays in parallel,
/// then merges them and does a final sort for a global suffix array.
fn build_parallel_suffix_array(genome: &str, num_chunks: usize) -> Vec<usize> {
    let chunk_size = (genome.len().max(1)) / num_chunks;
    // Create chunks as (start_index, chunk_string)
    let chunks: Vec<(usize, &str)> = (0..num_chunks)
        .map(|i| {
            let start = i * chunk_size;
            let end = if i == num_chunks - 1 {
                genome.len()
            } else {
                (i + 1) * chunk_size
            };
            (start, &genome[start..end])
        })
        .collect();

    // Build partial suffix arrays in parallel using Rayon
    let mut partial_results: Vec<(usize, Vec<usize>)> = chunks
        .par_iter()
        .map(|(offset, chunk)| {
            let sa = build_suffix_array(chunk);
            // Adjust each index by offset so they refer to the correct position in the full genome
            let adjusted_sa: Vec<usize> = sa.into_iter().map(|i| i + offset).collect();
            (*offset, adjusted_sa)
        })
        .collect();

    // Merge partial arrays into one vector of suffix starting points
    let mut combined: Vec<usize> = partial_results
        .iter_mut()
        .flat_map(|(_, sa)| sa.drain(..))
        .collect();

    // Finally, do a global sort to ensure suffixes are in correct lex order across chunk boundaries
    combined.sort_by(|&a, &b| genome[a..].cmp(&genome[b..]));
    combined
}

fn main() -> Result<()> {
    // 1) Load a large genomic sequence (in real HPC usage, might be GB-scale FASTA).
    let genome = load_genome("example_genome.fasta")?;

    // 2) Build a suffix array in parallel, splitting the sequence into multiple chunks.
    let num_chunks = 8; // depending on your HPC node or multi-core system
    let suffix_array = build_parallel_suffix_array(&genome, num_chunks);

    // 3) Display some stats.
    println!("Genome length: {}", genome.len());
    println!("Suffix array length: {}", suffix_array.len());
    println!("First 10 entries in suffix array: {:?}", &suffix_array[..10.min(suffix_array.len())]);

    Ok(())
}
