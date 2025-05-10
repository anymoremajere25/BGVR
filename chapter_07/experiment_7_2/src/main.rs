use clap::Parser;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::error::Error;
use std::fs::{self, File};
use std::io::{BufReader, BufWriter};

/// CLI arguments
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input reference FASTA file
    #[arg(long)]
    input: String,

    /// Output JSON file
    #[arg(long)]
    output: Option<String>,
}

/// Partial k-mer index
#[derive(Debug, Clone, Serialize, Deserialize)]
struct PartialIndex {
    chunk_id: usize,
    kmer_map: HashMap<String, usize>,
}

#[derive(Debug, Clone)]
struct IndexConfig {
    k: usize,
    chunk_size: usize,
    min_frequency: usize,
    filter_invalid: bool,
}

impl Default for IndexConfig {
    fn default() -> Self {
        Self {
            k: 31,
            chunk_size: 1_000_000,
            min_frequency: 1,
            filter_invalid: true,
        }
    }
}

fn is_valid_kmer(kmer: &str) -> bool {
    kmer.chars().all(|c| matches!(c, 'A' | 'C' | 'G' | 'T'))
}

fn build_partial_index(sequence: &str, chunk_id: usize, k: usize, filter: bool) -> PartialIndex {
    let mut kmer_map = HashMap::new();
    for i in 0..=sequence.len().saturating_sub(k) {
        let kmer = &sequence[i..i + k];
        if !filter || is_valid_kmer(kmer) {
            *kmer_map.entry(kmer.to_string()).or_insert(0) += 1;
        }
    }
    PartialIndex { chunk_id, kmer_map }
}

fn create_overlapping_chunks(sequence: &str, chunk_size: usize, k: usize) -> Vec<String> {
    if sequence.is_empty() || chunk_size == 0 {
        return Vec::new();
    }

    let overlap = k.saturating_sub(1);
    let mut chunks = Vec::new();
    let mut start = 0;

    while start < sequence.len() {
        let end = std::cmp::min(start + chunk_size, sequence.len());
        chunks.push(sequence[start..end].to_string());
        if end >= sequence.len() {
            break;
        }
        start = end.saturating_sub(overlap);
    }

    chunks
}

fn merge_partial_indexes(indexes: &[PartialIndex], min_freq: usize) -> HashMap<String, usize> {
    let mut global_map = HashMap::new();
    for idx in indexes {
        for (kmer, count) in &idx.kmer_map {
            *global_map.entry(kmer.clone()).or_insert(0) += count;
        }
    }
    if min_freq > 1 {
        global_map.retain(|_, count| *count >= min_freq);
    }
    global_map
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().collect();

    if args.contains(&"--merge".to_string()) {
        let merge_index = args.iter().position(|x| x == "--merge").unwrap();
        let output_index = args.iter().position(|x| x == "--output").unwrap();

        let input_files = &args[merge_index + 1..output_index];
        let output_file = &args[output_index + 1];

        let mut merged_map = HashMap::new();
        for path in input_files {
            let file = File::open(path)?;
            let reader = BufReader::new(file);
            let partial: HashMap<String, usize> = serde_json::from_reader(reader)?;
            for (k, v) in partial {
                *merged_map.entry(k).or_insert(0) += v;
            }
        }

        let writer = BufWriter::new(File::create(output_file)?);
        serde_json::to_writer(writer, &merged_map)?;
        println!("Merged index written to {}", output_file);
        return Ok(());
    }

    // normal single-file mode
    let input_index = args.iter().position(|x| x == "--input").unwrap();
    let output_index = args.iter().position(|x| x == "--output").unwrap();

    let input_path = &args[input_index + 1];
    let output_path = &args[output_index + 1];

    let config = IndexConfig::default();
    let reference = fs::read_to_string(input_path)?;

    let chunks = create_overlapping_chunks(&reference, config.chunk_size, config.k);
    let partial_indexes: Vec<_> = chunks
        .par_iter()
        .enumerate()
        .map(|(i, chunk)| build_partial_index(chunk, i, config.k, config.filter_invalid))
        .collect();

    let global_index = merge_partial_indexes(&partial_indexes, config.min_frequency);

    let file = File::create(output_path)?;
    let writer = BufWriter::new(file);
    serde_json::to_writer(writer, &global_index)?;

    println!("Global index written to {}", output_path);
    Ok(())
}