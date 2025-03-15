## 4.2 Graph-Based Models for Gene Regulatory Networks (GRNs)

project_42_1

### Explanation of the Rust Program: Computing a Correlation-Based Gene Network

This Rust program efficiently computes a correlation-based adjacency matrix for gene expression data using parallel computing. The program generates synthetic gene expression data, calculates Pearson correlation coefficients between genes, and stores the resulting adjacency matrix in a binary file.
### 1. Input Parameters and Execution Flow

The program accepts command-line arguments to configure:

    --num-genes → Number of genes (default: 1000)

    --num-samples → Number of samples (default: 50)

    --output → Name of the output file (default: partial_adjacency.bin)

If no arguments are provided, it runs with the default settings.
Execution Steps:

    Parse command-line arguments to get the number of genes, number of samples, and output filename.

    Generate synthetic gene expression data (a random matrix of gene expression values).

    Compute Pearson correlation coefficients for all gene pairs.

    Store the adjacency matrix in a binary file.

### 2. Generating Synthetic Gene Expression Data

The function generate_synthetic_expression(num_genes, num_samples) creates a 2D array (Array2<f64>) of shape (num_genes, num_samples). Each value in this matrix represents the expression level of a gene in a specific sample, randomly drawn from the range 0.0 to 1000.0.

fn generate_synthetic_expression(num_genes: usize, num_samples: usize) -> Array2<f64> {
    let mut data = Array2::<f64>::zeros((num_genes, num_samples));
    let mut rng = thread_rng(); // Random number generator

    for i in 0..num_genes {
        for j in 0..num_samples {
            data[[i, j]] = rng.gen_range(0.0..1000.0); // Random values between 0 and 1000
        }
    }
    data
}

Purpose: Simulates real-world gene expression data, which typically varies across samples.
###  3. Computing Pearson Correlation Between Genes

The function pearson_correlation(x, y) calculates the correlation between two gene expression profiles (two rows from the matrix).
Pearson Correlation Formula:
r=∑(xi−xˉ)(yi−yˉ)∑(xi−xˉ)2∑(yi−yˉ)2
r=∑(xi​−xˉ)2
​∑(yi​−yˉ​)2
​∑(xi​−xˉ)(yi​−yˉ​)​

    x and y → Expression values of two genes across samples.

    \bar{x} and \bar{y} → Mean expression levels of genes.

fn pearson_correlation(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len();
    let mean_x = x.iter().sum::<f64>() / n as f64;
    let mean_y = y.iter().sum::<f64>() / n as f64;

    let mut num = 0.0;
    let mut den_x = 0.0;
    let mut den_y = 0.0;
    for i in 0..n {
        let dx = x[i] - mean_x;
        let dy = y[i] - mean_y;
        num += dx * dy;
        den_x += dx * dx;
        den_y += dy * dy;
    }
    let denom = den_x.sqrt() * den_y.sqrt();
    if denom == 0.0 { 0.0 } else { num / denom }
}

Purpose: Measures how strongly the expression patterns of two genes are related.
#### 4. Parallel Computation of the Adjacency Matrix

Since calculating correlations for 1000 genes (~500,000 pairs) is computationally expensive, Rayon is used for parallel processing.
Steps:

    A shared adjacency matrix (Array2<f64>) is created with all values initialized to zero.

    Parallel iteration over gene pairs using par_iter().

    Each (i, j) pair is computed only once (upper triangle of the matrix).

    Mutex (Arc<Mutex<Array2<f64>>) ensures safe access to the shared adjacency matrix.

(0..num_genes).into_par_iter().for_each(|i| {
    let row_i = expression_data.slice(s![i, ..]);
    for j in (i+1)..num_genes {
        let row_j = expression_data.slice(s![j, ..]);
        let corr = pearson_correlation(&row_i.to_vec(), &row_j.to_vec());

        let mut adj_mut = adjacency.lock().unwrap();
        adj_mut[[i, j]] = corr;
        adj_mut[[j, i]] = corr;
    }
});

Purpose: Computes correlation efficiently by processing different gene pairs in parallel.
### 5. Writing the Correlation Matrix to a Binary File

After computation, the adjacency matrix is written to a binary file (partial_adjacency.bin). Each value is stored as native-endian floating-point bytes.

let file = File::create(&output_file)?;
let mut writer = BufWriter::new(file);
for i in 0..num_genes {
    for j in 0..num_genes {
        writer.write_all(&final_adj_matrix[[i, j]].to_ne_bytes())?;
    }
}
println!("Correlation adjacency matrix written to {}", output_file);

Purpose: Saves space and allows fast loading of the matrix in future computations.
### 6. Example Output

If the program runs correctly with default settings:

Number of genes: 1000
Number of samples: 50
Output file: partial_adjacency.bin
Correlation adjacency matrix written to partial_adjacency.bin

Output File (partial_adjacency.bin):

    Contains the correlation matrix in binary format.

    Can be used in further gene network analysis.

#### 7. Summary
Step	Description
Generate Synthetic Data	Create a matrix of random gene expression values.
Compute Pearson Correlation	Measure relationships between gene expression patterns.
Use Parallel Processing	Speed up computation with Rayon.
Store Adjacency Matrix	Save results in a binary file for future analysis.
Key Features:

✅ Efficient Parallel Computation using Rayon
✅ Realistic Synthetic Data with random gene expression levels
✅ Binary File Output for easy storage and retrieval
