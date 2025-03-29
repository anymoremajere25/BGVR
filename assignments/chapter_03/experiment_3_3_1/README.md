### 3.3. Graph-Based Data Structures for Genome Assembly and Beyond  
**experiment_3_3_1**  

The following Rust code demonstrates the construction of partial de Bruijn graphs from FASTQ data. It utilizes **needletail** for efficient sequence parsing and **Petgraph** to model overlapping k-mers as an undirected graph. By splitting sequencing reads into manageable chunks and processing them in parallel with **Rayon**, each segment contributes to a partial de Bruijn graph. These partial graphs are then serialized for later merging or further analysis. This modular approach enables efficient handling of large genomic datasets by leveraging parallel computing and optimizing memory usage.  

### Implementation Overview  

The program reads sequences from a FASTQ file, partitions them into equal-sized chunks, and assigns each chunk to a separate thread. Within each chunk, k-mers are extracted from sequencing reads, added as nodes in the graph, and linked by edges representing consecutive overlaps. After constructing the local subgraph, the data (nodes and edges) is serialized into JSON format for modular processing. This approach allows intermediate results to be stored and combined later into a global de Bruijn graph.  

**Project Directory Structure:**  
```
experiment_3_3_1/
    Cargo.toml  (Dependency configuration)
    experiment_3_3_1/src/
        main.rs  (Rust implementation)
        reads.fq.rar  (Compressed sequencing reads)
        partial_debrujin_graphs.json.rar  (Compressed JSON output)
        output.txt  (Program output)
```
### Running the Program  

Execute the following command in PowerShell:  
```sh
cargo run | tee output.txt
```
This runs `main.rs` and saves the output in `output.txt`.  

### Dependencies  
```toml
[dependencies]
rayon = "1.7"
needletail = "0.6"
serde = { version = "1", features = ["derive"] }
serde_json = "1"
petgraph = "0.7.1"
```
---

## **Understanding the Output**  

### **1. `output.txt`**  
```
Wrote partial de Bruijn graphs to partial_debruijn_graphs.json
```
This message confirms that the de Bruijn graph construction was completed successfully, and the results were saved in `partial_debruijn_graphs.json`.  

### **2. `partial_debruijn_graphs.json`**  
This JSON file contains the generated partial de Bruijn graphs with the following structure:  

- **nodes**: List of unique k-mers (31-mers).  
- **edges**: List of tuples `(i, j)`, where an edge exists between `nodes[i]` and `nodes[j]`.  
- **k**: The k-mer length, set to **31**.  

**Example JSON Output:**  
```json
{
  "nodes": [
    "ACGGAGGATGCGAGCGTTATCCGGATTTATT",
    "CGGAGGATGCGAGCGTTATCCGGATTTATTG",
    "GGAGGATGCGAGCGTTATCCGGATTTATTGG",
    "GAGGATGCGAGCGTTATCCGGATTTATTGGG",
    "AGGATGCGAGCGTTATCCGGATTTATTGGGT",
    "GGATGCGAGCGTTATCCGGATTTATTGGGTT"
  ],
  "edges": [
    [0, 1],
    [1, 2],
    [2, 3],
    [3, 4],
    [4, 5]
  ],
  "k": 31
}
```
### **Understanding the Graph Representation**  

- Each **node** represents a unique k-mer of length **31**.  
- Each **edge** indicates overlap between consecutive k-mers (i.e., sharing `k-1` characters).  
- This structure enables genome assembly by linking overlapping k-mers.  

---

## **Algorithm Breakdown**  

### **Step 1: Read Input FASTQ/FASTA File**  
- The program reads sequencing data from `reads.fq`.  
- Uses **needletail** for efficient parsing of FASTA/FASTQ formats.  

### **Step 2: Partition Reads into Chunks**  
- Reads are split into **chunks of 100,000 sequences** each.  
- Enables parallel processing for better performance.  

### **Step 3: Construct Partial de Bruijn Graphs**  
For each chunk, the program:  
1. Extracts **31-mers** from reads.  
2. Adds each **k-mer** as a node in the graph.  
3. Creates edges between **consecutive k-mers**.  
4. Stores the partial graph.  

### **Step 4: Parallel Processing with Rayon**  
- Uses **Rayon** to process chunks concurrently.  
- Each thread constructs a **partial de Bruijn graph** independently.  

### **Step 5: Output Partial Graphs**  
- Each partial graph is **serialized to JSON (`partial_debruijn_graphs.json`)**.  
- Allows further processing, such as merging into a complete **global de Bruijn graph**.  

---

## **Key Takeaways**  

### ✅ **Successful Graph Construction**  
- The program correctly generates partial de Bruijn graphs, as confirmed by the JSON output.  

### ✅ **Efficient Parallelization**  
- The use of **Rayon** speeds up processing by parallelizing graph construction.  

### ✅ **Scalability**  
- A **chunk-based approach** prevents memory overload, making it feasible to process large datasets.  

### ✅ **Correct Graph Representation**  
- The edges accurately represent **overlapping k-mers**, ensuring proper connectivity.  

### **Next Steps:**  
- Merge partial graphs into a **complete de Bruijn graph**.  
- Apply the method for **genome assembly, variant detection, or error correction**.  

### **Final Thoughts**  
This Rust-based pipeline efficiently constructs **partial de Bruijn graphs** from sequencing data.  
- **Parallel execution** optimizes performance.  
- **Graph-based modeling** supports genome assembly and sequence analysis.  
- The approach is modular, scalable, and well-suited for **large genomic datasets**.
