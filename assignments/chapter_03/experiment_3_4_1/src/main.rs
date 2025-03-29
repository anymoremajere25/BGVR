use rayon::prelude::*;
use needletail::parse_fastx_file;
use serde::{Serialize, Deserialize};
use std::path::PathBuf;
use std::fs;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::error::Error;

/// A Bloom filter for k-mers, supporting partial merges in HPC settings.
#[derive(Serialize, Deserialize)]
struct BloomFilter {
    bits: Vec<u8>,
    num_bits: usize,
    num_hashes: usize,
    k: usize,
}

impl BloomFilter {
    pub fn new(num_bits: usize, num_hashes: usize, k: usize) -> Self {
        let byte_len = (num_bits + 7) / 8;
        BloomFilter {
            bits: vec![0; byte_len],
            num_bits,
            num_hashes,
            k,
        }
    }

    /// Simple FNV-based hash, seeded by `seed`.
    fn hash_kmer(kmer: &str, seed: u64) -> u64 {
        let mut hash_val = 0xcbf29ce484222325 ^ seed;
        for b in kmer.as_bytes() {
            hash_val ^= *b as u64;
            hash_val = hash_val.wrapping_mul(0x100000001b3);
        }
        hash_val
    }

    /// Insert a k-mer into the Bloom filter.
    pub fn insert(&mut self, kmer: &str) {
        for i in 0..self.num_hashes {
            let h = Self::hash_kmer(kmer, i as u64);
            let bit_index = (h as usize) % self.num_bits;
            let byte_idx = bit_index / 8;
            let bit_idx = bit_index % 8;
            self.bits[byte_idx] |= 1 << bit_idx;
        }
    }

    /// Check membership of a k-mer.
    pub fn contains(&self, kmer: &str) -> bool {
        for i in 0..self.num_hashes {
            let h = Self::hash_kmer(kmer, i as u64);
            let bit_index = (h as usize) % self.num_bits;
            let byte_idx = bit_index / 8;
            let bit_idx = bit_index % 8;
            if (self.bits[byte_idx] & (1 << bit_idx)) == 0 {
                return false;
            }
        }
        true
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    // 1) Parse command-line arguments
    let input_path = std::env::args().nth(1).unwrap_or_else(|| "reads.fq".into());
    let k: usize = std::env::args().nth(2).unwrap_or_else(|| "31".into()).parse()?;

    // 2) Initialize Bloom filter
    let num_bits: usize = 10_000_000;
    let num_hashes: usize = 3;
    let mut bloom = BloomFilter::new(num_bits, num_hashes, k);

    // 3) Read FASTQ/FASTA input
    let mut reader = parse_fastx_file(PathBuf::from(&input_path))?;
    let mut all_reads = Vec::new();
    while let Some(record) = reader.next() {
        let seqrec = record?;
        all_reads.push(seqrec.seq().to_vec());
    }

    // 4) Gather all k-mers in parallel (no mutation of bloom)
    let total_kmers = AtomicUsize::new(0);

    let all_kmers: Vec<String> = all_reads
        .par_iter()
        .flat_map_iter(|read| {
            let mut local_kmers = Vec::new();
            for i in 0..=read.len().saturating_sub(k) {
                let kmer_str = String::from_utf8_lossy(&read[i..i + k]).to_string();
                local_kmers.push(kmer_str);
            }
            total_kmers.fetch_add(local_kmers.len(), Ordering::Relaxed);
            local_kmers
        })
        .collect();

    // 5) Insert k-mers into Bloom filter sequentially
    for kmer in &all_kmers {
        bloom.insert(kmer);
    }

    // Demonstrate usage of the `contains` method
    let test_kmer = if let Some(kmer) = all_kmers.get(0) {
        kmer
    } else {
        "ATG" // fallback if no reads
    };
    println!(
        "Bloom contains first k-mer '{}': {}",
        test_kmer,
        bloom.contains(test_kmer)
    );

    // 6) Serialize partial Bloom filter
    let serialized = serde_json::to_string_pretty(&bloom)?;
    fs::write("bloom.json", serialized)?;

    println!(
        "Constructed Bloom filter ({} k-mers processed), result written to bloom.json",
        total_kmers.load(Ordering::Relaxed)
    );
    Ok(())
}
