use std::collections::HashMap;
use rayon::prelude::*;
use bio::io::fasta;

/// Represents a node in the pangenome graph, storing a k-mer label.
#[derive(Clone, Debug, Eq, PartialEq, Hash)]
struct PGNode {
    kmer: String,
}

/// An adjacency list mapping each node to its successor nodes.
type PGGraph = HashMap<PGNode, Vec<PGNode>>;

/// Builds a pangenome graph by merging k-mers from multiple FASTA files.
/// Each file corresponds to one haplotype/genome. Overlapping k-mers
/// are treated as edges, and equivalent k-mers across files are merged.
fn build_pangenome_graph(k: usize, fasta_files: &[&str]) -> PGGraph {
    let partial_graphs: Vec<PGGraph> = fasta_files
        .par_iter()
        .map(|path| {
            println!("Reading file: {}", path); // Debug: indicate which file is being processed
            let reader = fasta::Reader::from_file(path)
                .unwrap_or_else(|_| panic!("Cannot open FASTA file: {}", path));
            let mut local_graph = PGGraph::new();
            for record in reader.records() {
                let rec = record.expect("Invalid FASTA record");
                let seq_bytes = rec.seq();
                let seq_str = String::from_utf8(seq_bytes.to_vec())
                    .expect("Non-UTF8 sequence data");

                println!("Processing sequence: {}", seq_str); // Debug: print sequence being processed

                // If the sequence is shorter than k, skip it.
                if seq_str.len() < k {
                    println!("Skipping sequence, length < k: {}", seq_str); // Debug: skipped due to length
                    continue;
                }

                // Build the k-mer graph for the sequence.
                for i in 0..seq_str.len().saturating_sub(k) {
                    let node_kmer = &seq_str[i..i + k];
                    let edge_kmer = &seq_str[i + 1..i + k + 1];

                    let node = PGNode {
                        kmer: node_kmer.to_owned(),
                    };
                    let edge_node = PGNode {
                        kmer: edge_kmer.to_owned(),
                    };

                    local_graph
                        .entry(node)
                        .or_insert_with(Vec::new)
                        .push(edge_node);
                }
            }
            local_graph
        })
        .collect();

    // Merge all partial graphs into a single pangenome graph.
    partial_graphs
        .into_par_iter()
        .reduce(
            || PGGraph::new(),
            |mut acc, local| {
                for (node, successors) in local {
                    acc.entry(node)
                        .or_insert_with(Vec::new)
                        .extend(successors);
                }
                acc
            },
        )
}

fn main() {
    let haplotypes = &["src/haplotype1.fasta", "src/haplotype2.fasta", "src/haplotype3.fasta"];
    let k = 10;
    let pangenome_graph = build_pangenome_graph(k, haplotypes);

    println!(
        "Constructed a pangenome graph with {} nodes.",
        pangenome_graph.len()
    );

    // Print a small subset of the graph
    for (node, edges) in pangenome_graph.iter().take(5) {
        println!("Node: {} -> {:?}", node.kmer, edges.iter().map(|e| &e.kmer).collect::<Vec<_>>());
    }
}
