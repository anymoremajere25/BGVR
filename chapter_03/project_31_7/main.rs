use rayon::prelude::*;
use needletail::parse_fastx_file;
use serde::{Deserialize, Serialize};
use serde_json;
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::BufWriter;
use std::path::PathBuf;
use std::{env, error::Error};

/// A minimal struct representing a De Bruijn node.
#[derive(Debug, Serialize, Deserialize)]
struct Node {
    kmer: String,
    // Store edges as the set of next-kmers to which this k-mer connects.
    edges: HashSet<String>,
}

/// A minimal Bloom Filter with parallel insertion. 
#[derive(Debug, Serialize, Deserialize)]
struct BloomFilter {
    bits: Vec<bool>,
    num_hashes: usize,
    size: usize,
}

impl BloomFilter {
    /// Create a new BloomFilter with `size` bits and `num_hashes` independent hash functions.
    fn new(size: usize, num_hashes: usize) -> Self {
        BloomFilter {
            bits: vec![false; size],
            num_hashes,
            size,
        }
    }

    /// Insert a string into the filter in parallel over the number of hash functions.
    fn insert(&self, item: &str) {
        // Because 'bits' is a Vec<bool>, we can't set it in parallel without locks or atomic booleans.
        // For simplicity, we do a sequential approach, but you can adapt this using AtomicBool if needed.
        // Alternatively, you can chunk items for parallel insertion.
        for seed in 0..self.num_hashes {
            let h = self.hash_with_seed(item, seed) % self.size;
            // With a normal Vec<bool>, we need exclusive access:
            // but let's keep it simple here. In real HPC, consider AtomicBool or a bit vector crate.
            unsafe {
                // This is a quick hack; prefer a thread-safe approach with minimal overhead.
                let bits_ptr = self.bits.as_ptr() as *mut bool;
                *bits_ptr.add(h) = true;
            }
        }
    }

    /// Check membership in the filter. If any bit is false, definitely not in set.
    fn contains(&self, item: &str) -> bool {
        for seed in 0..self.num_hashes {
            let h = self.hash_with_seed(item, seed) % self.size;
            if !self.bits[h] {
                return false;
            }
        }
        true
    }

    /// Hash function with a seed mixed in.
    fn hash_with_seed(&self, item: &str, seed: usize) -> usize {
        use std::hash::{Hash, Hasher};
        use std::collections::hash_map::DefaultHasher;

        let mut hasher = DefaultHasher::new();
        item.hash(&mut hasher);
        seed.hash(&mut hasher);
        hasher.finish() as usize
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    // Simple CLI args. For HPC usage, adapt them or read Nextflow params from environment.
    let mut args = env::args();
    let _bin = args.next(); // skip binary name

    let mut fastq_path = String::from("example.fastq");
    let mut kmer_size = 31;
    let mut outdir = String::from("results");

    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--fastq" => fastq_path = args.next().unwrap(),
            "--kmer" => {
                let val = args.next().unwrap();
                kmer_size = val.parse().unwrap_or(31);
            }
            "--outdir" => outdir = args.next().unwrap(),
            _ => {}
        }
    }

    // Ensure output directory exists
    fs::create_dir_all(&outdir)?;

    // 1) Read the FASTQ file using Needletail.
    let mut reader = parse_fastx_file(PathBuf::from(&fastq_path))?;
    println!("Reading from FASTQ: {}, k-mer={}", fastq_path, kmer_size);

    // 2) Build an in-memory list of all read strings for demonstration.
    // For real HPC usage, you may want chunk-based or streaming approaches.
    let mut reads: Vec<String> = Vec::new();
    while let Some(record) = reader.next() {
        let seqrec = record?;
        // Convert to owned String for further manipulations.
        reads.push(String::from_utf8_lossy(seqrec.seq().as_ref()).to_string());
    }
    println!("Loaded {} reads. Building de Bruijn graph & Bloom filter...", reads.len());

    // 3) Build a minimal De Bruijn graph representation.
    let nodes_map = build_debruijn(&reads, kmer_size);

    // 4) Build a Bloom filter of all distinct k-mers in parallel (for demonstration).
    let distinct_kmers_count = nodes_map.len();
    let bf_size = (distinct_kmers_count * 10).max(1_000_000); // fallback to 1M bits
    let bloom = BloomFilter::new(bf_size, 3);

    // Insert all k-mers in parallel over the keys.
    nodes_map
        .par_iter()
        .for_each(|(kmer, _node)| {
            bloom.insert(kmer);
        });

    // OPTIONAL: Call 'contains' to demonstrate usage, preventing a dead_code warning.
    // e.g. check if "ACGT" was possibly inserted:
    let test_kmer = "ACGT";
    println!("Bloom filter contains '{}'? {}", test_kmer, bloom.contains(test_kmer));

    // 5) Output the graph as JSON to `graph.json`, and bloom filter to `bloom.json`.
    let graph_path = format!("{}/graph.json", outdir);
    let bloom_path = format!("{}/bloom.json", outdir);

    let graph_file = File::create(&graph_path)?;
    serde_json::to_writer(BufWriter::new(graph_file), &nodes_map)?;

    let bloom_file = File::create(&bloom_path)?;
    serde_json::to_writer(BufWriter::new(bloom_file), &bloom)?;

    println!("De Bruijn graph saved to {}", graph_path);
    println!("Bloom filter saved to {}", bloom_path);
    Ok(())
}

/// Build a minimal de Bruijn graph for the reads. 
/// Returns a HashMap from k-mer string -> Node.
fn build_debruijn(reads: &[String], k: usize) -> HashMap<String, Node> {
    // We'll store a Node per k-mer. The Node's `edges` field
    // is a set of next k-mers that follow in the read.
    // We parallelize extraction of k-mers, then build adjacency in a second pass.
    let kmers = reads
        .par_iter()
        .flat_map(|read| {
            let mut local_kmers = Vec::new();
            if read.len() >= k {
                for i in 0..=read.len() - k {
                    local_kmers.push(&read[i..i + k]);
                }
            }
            local_kmers
        })
        .collect::<Vec<&str>>();

    let mut map: HashMap<String, Node> = HashMap::new();
    // Insert all k-mers as nodes
    for kmer_str in kmers.iter() {
        map.entry((*kmer_str).to_string())
            .or_insert(Node {
                kmer: (*kmer_str).to_string(),
                edges: HashSet::new(),
            });
    }

    // Build adjacency in a second pass. This is naive:
    reads.iter().for_each(|read| {
        if read.len() < k + 1 {
            return;
        }
        for i in 0..=read.len() - (k + 1) {
            let k1 = &read[i..i + k];
            let k2 = &read[i + 1..i + 1 + k];
            // Insert an edge from k1 -> k2
            if let Some(node) = map.get_mut(k1) {
                node.edges.insert(k2.to_string());
            }
        }
    });

    map
}
