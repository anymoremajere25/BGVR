## 4.6 SINGLE-CELL FUNCTIONAL GENOMICS
### Explanation of the k-NN Graph Construction and Pseudotime Computation [project_46_1]

This Rust program builds a k-nearest neighbor (k-NN) graph from single-cell expression data and estimates pseudotime using Breadth-First Search (BFS). Let's go step by step through the process:
#### 1. Data Representation

We define two main data structures:

    Cell – Represents a single biological cell with an expression profile (a vector of floating-point numbers).

struct Cell {
    _id: usize, // Unique identifier
    expression: Vec<f64>, // Gene expression values
}

KnnGraph – A graph where each node represents a cell, and edges connect it to its k-nearest neighbors.

    struct KnnGraph {
        edges: Vec<Vec<usize>>, // Adjacency list representation
    }

#### 2. Constructing the k-NN Graph
#### Step 1: Compute Pairwise Distances

    Each cell’s expression vector is compared with all other cells to compute distances.
    The Euclidean distance formula is used:
    d(a,b)=∑(ai−bi)2
    d(a,b)=∑(ai​−bi​)2

​

    fn euclidean_dist(a: &[f64], b: &[f64]) -> f64 {
        a.iter().zip(b.iter()).map(|(x, y)| (x - y).powi(2)).sum::<f64>().sqrt()
    }

#### Step 2: Select k Nearest Neighbors

    Each cell finds its k closest cells (excluding itself).
    This is done in parallel using Rayon, a Rust library for parallel computing.

    let edges_arc = Arc::new(Mutex::new(vec![Vec::new(); n])); // Shared memory for graph edges
    (0..n).into_par_iter().for_each(|i| {
        let mut dists = Vec::new();
        for j in 0..n {
            if i != j {
                let dist = euclidean_dist(&cells[i].expression, &cells[j].expression);
                dists.push((j, dist)); // Store index and distance
            }
        }
        dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap()); // Sort by distance
        let top_k: Vec<usize> = dists.iter().take(k).map(|&(idx, _)| idx).collect();
        edges_arc.lock().unwrap()[i] = top_k;
    });

    Parallelism helps process large datasets faster by distributing distance computations across CPU cores.

#### 3. Computing Pseudotime

Once the k-NN graph is built, we estimate pseudotime. Pseudotime is a measure of how far each cell is from a chosen root cell, calculated as the shortest distance using BFS.
#### Step 1: Initialize BFS

    We start from a root cell (cell 0).
    pseudotime values are initialized to ∞ (unvisited).
    The root cell is set to 0.0 (starting point).

    let mut pseudotime = vec![f64::INFINITY; n];
    pseudotime[root] = 0.0;
    let mut queue = VecDeque::new();
    queue.push_back(root);

#### Step 2: Traverse the Graph with BFS

    BFS explores the nearest neighbors first, assigning each an incremented time value.
    For each cell:
        Visit all connected neighbors.
        If a neighbor is unvisited (still has ∞ pseudotime), assign it a value current_time + 1.
        Push the neighbor to the queue for further processing.

    while let Some(current) = queue.pop_front() {
        let current_time = pseudotime[current];
        for &nbr in &graph.edges[current] {
            if pseudotime[nbr].is_infinite() {
                pseudotime[nbr] = current_time + 1.0;
                queue.push_back(nbr);
            }
        }
    }

    This algorithm ensures that closer cells get smaller pseudotime values, simulating a temporal progression.

#### 4. Output Interpretation
Example k-NN Graph

k-NN edges:
Cell 0 neighbors = [4, 1]
Cell 1 neighbors = [4, 0]
Cell 2 neighbors = [3, 1]
Cell 3 neighbors = [2, 1]
Cell 4 neighbors = [0, 1]

    Each cell is connected to its two nearest neighbors (k=2).

Pseudotime Output

Pseudotime from root=0: [0.0, 1.0, inf, inf, 1.0]

    Cells 2 and 3 have inf pseudotime → They are disconnected from Cell 0!
    This means the k-NN graph is not fully connected, so BFS cannot reach every node.

### 5. Summary

✔ Parallelized k-NN Construction → Efficient for large datasets
✔ Graph-Based Pseudotime Estimation → BFS-based approach
⚠ Disconnected Graph Issue → If k is too small, some cells might not connect to the main network
#### Possible Improvements

    Increase k → Helps prevent disconnected components.
    Check Connectivity → Ensure all cells are reachable before running BFS.
