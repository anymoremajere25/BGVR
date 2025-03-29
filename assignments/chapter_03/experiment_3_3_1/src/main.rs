use rayon::prelude::*;
use needletail::parse_fastx_file;
use petgraph::graph::Graph;
use petgraph::Undirected;
use serde::{Serialize, Deserialize};
use std::path::PathBuf;
use std::fs;
use std::collections::HashMap;
use std::error::Error;

/// This struct captures a subgraph of the entire de Bruijn representation.
#[derive(Debug, Serialize, Deserialize)]
struct PartialDeBruijnGraph {
    nodes: Vec<String>,
    edges: Vec<(usize, usize)>,
    k: usize,
}

fn main() -> Result<(), Box<dyn Error>> {
    let reads_path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "reads.fq".to_string());
    let k: usize = std::env::args()
        .nth(2)
        .unwrap_or_else(|| "31".to_string())
        .parse()?;
    let chunk_size = 100_000;

    // Collect sequences from the FASTQ/FASTA file.
    let mut reads_data = Vec::new();
    let mut reader = parse_fastx_file(PathBuf::from(&reads_path))?;
    while let Some(record) = reader.next() {
        let seqrec = record?;
        // Convert from Cow<'_, [u8]> to a String.
        let seq_str = String::from_utf8_lossy(seqrec.seq().as_ref()).to_string();
        reads_data.push(seq_str);
    }

    // Partition reads into chunks for parallel building of partial graphs.
    let reads_chunks: Vec<_> = reads_data
        .chunks(chunk_size)
        .map(|c| c.to_vec())
        .collect();

    // Build partial De Bruijn graphs concurrently.
    let partial_graphs: Vec<PartialDeBruijnGraph> = reads_chunks
        .par_iter()
        .map(|chunk_reads| build_partial_debruijn(&chunk_reads, k))
        .collect();

    // Serialize all partial graphs to JSON and write to disk.
    let serialized = serde_json::to_string_pretty(&partial_graphs)?;
    fs::write("partial_debruijn_graphs.json", serialized)?;
    println!("Wrote partial de Bruijn graphs to partial_debruijn_graphs.json");

    Ok(())
}

/// Constructs a subgraph (partial De Bruijn) from a collection of reads.
fn build_partial_debruijn(reads: &Vec<String>, k: usize) -> PartialDeBruijnGraph {
    // Use new_undirected() to build an undirected graph when specifying the Undirected type.
    let mut graph = Graph::<String, (), Undirected>::new_undirected();
    let mut node_map = HashMap::new();

    for read in reads {
        if read.len() < k {
            continue;
        }
        // For each k-mer in the read...
        for i in 0..=(read.len() - k) {
            let kmer = &read[i..i + k];
            // Insert node if not present
            let node_index = *node_map.entry(kmer.to_string()).or_insert_with(|| {
                graph.add_node(kmer.to_string())
            });
            // Link to the previous k-mer, forming edges between consecutive k-mers.
            if i > 0 {
                let prev_kmer = &read[i - 1..i - 1 + k];
                if let Some(&prev_index) = node_map.get(prev_kmer) {
                    graph.add_edge(prev_index, node_index, ());
                }
            }
        }
    }

    // Collect node labels and edges for serialization.
    let mut nodes = Vec::new();
    let mut edges = Vec::new();

    // Extract node labels (k-mers) in their internal order.
    for idx in graph.node_indices() {
        nodes.push(graph[idx].clone());
    }

    // Map each NodeIndex to a position in the nodes vector.
    let mut node_idx_map = HashMap::new();
    for (i, idx) in graph.node_indices().enumerate() {
        node_idx_map.insert(idx, i);
    }

    // Collect edges as pairs of node indices (s_idx, t_idx).
    for edge_idx in graph.edge_indices() {
        if let Some((src, tgt)) = graph.edge_endpoints(edge_idx) {
            let s_idx = node_idx_map[&src];
            let t_idx = node_idx_map[&tgt];
            edges.push((s_idx, t_idx));
        }
    }

    PartialDeBruijnGraph { nodes, edges, k }
}
