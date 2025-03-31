use rayon::prelude::*;
use std::collections::VecDeque; // We keep VecDeque for BFS
use std::sync::{Arc, Mutex};

/// A simple structure to hold single-cell data, e.g., expression vectors.
#[derive(Clone)]
struct Cell {
    // Renamed `id` to `_id` to avoid warnings if not used.
    _id: usize,
    expression: Vec<f64>,
}

/// Represents a k-nearest neighbor graph for cells. `edges[id]` contains the neighbors of cell `id`.
#[derive(Debug)]
struct KnnGraph {
    edges: Vec<Vec<usize>>,
}

/// Constructs a k-NN graph in parallel by computing distances between cells and picking the top `k`.
/// For HPC-scale datasets, consider approximate methods (e.g., Annoy, Faiss) or chunk-based merges.
fn build_knn_graph(cells: &[Cell], k: usize) -> KnnGraph {
    let n = cells.len();
    let edges_arc = Arc::new(Mutex::new(vec![Vec::new(); n]));

    // We compute distances for each cell -> all others in parallel, then pick top k.
    // For extremely large `n`, a local partial merges approach can reduce lock contention.
    (0..n).into_par_iter().for_each(|i| {
        let mut dists = Vec::with_capacity(n);
        for j in 0..n {
            if i == j {
                continue;
            }
            let dist = euclidean_dist(&cells[i].expression, &cells[j].expression);
            dists.push((j, dist));
        }
        // Sort by distance, then pick the top k.
        dists.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        let top_k: Vec<usize> = dists.iter().take(k).map(|&(idx, _)| idx).collect();

        let mut lock = edges_arc.lock().unwrap();
        lock[i] = top_k;
    });

    let final_edges = Arc::try_unwrap(edges_arc).unwrap().into_inner().unwrap();
    KnnGraph { edges: final_edges }
}

/// Example function to approximate "pseudotime" by BFS from a chosen root cell.
/// Each cell's pseudotime is the BFS distance from the root in the k-NN graph.
fn compute_pseudotime(graph: &KnnGraph, root: usize) -> Vec<f64> {
    let n = graph.edges.len();
    let mut pseudotime = vec![f64::INFINITY; n];
    let mut queue = VecDeque::new();

    pseudotime[root] = 0.0;
    queue.push_back(root);

    while let Some(current) = queue.pop_front() {
        let current_time = pseudotime[current];
        for &nbr in &graph.edges[current] {
            // If not visited (INFINITY), assign BFS distance = current_time + 1
            if pseudotime[nbr].is_infinite() {
                pseudotime[nbr] = current_time + 1.0;
                queue.push_back(nbr);
            }
        }
    }
    pseudotime
}

/// Basic Euclidean distance for demonstration.
fn euclidean_dist(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| (x - y).powi(2)).sum::<f64>().sqrt()
}

fn main() {
    // Suppose HPC ephemeral containers each handle part of a large single-cell dataset.
    // Here, we show a small synthetic set of cells for demonstration.
    let cells = vec![
        Cell {
            _id: 0,
            expression: vec![1.0, 2.0, 3.0],
        },
        Cell {
            _id: 1,
            expression: vec![2.0, 2.1, 3.2],
        },
        Cell {
            _id: 2,
            expression: vec![4.0, 5.0, 6.0],
        },
        Cell {
            _id: 3,
            expression: vec![3.9, 5.1, 6.1],
        },
        Cell {
            _id: 4,
            expression: vec![1.1, 2.2, 2.9],
        },
    ];

    // Build a k-NN graph in parallel. For HPC or cloud usage, partial merges
    // or approximate methods are recommended if `n` is extremely large.
    let k = 2;
    let graph = build_knn_graph(&cells, k);

    // Choose an arbitrary root cell (0) to define "pseudotime" as BFS distance from root.
    let pseudo_times = compute_pseudotime(&graph, 0);

    println!("k-NN edges:");
    for (i, nbrs) in graph.edges.iter().enumerate() {
        println!("Cell {} neighbors = {:?}", i, nbrs);
    }

    // Each cell's BFS distance from root can be considered a naive pseudotime
    println!("Pseudotime from root=0: {:?}", pseudo_times);
}
