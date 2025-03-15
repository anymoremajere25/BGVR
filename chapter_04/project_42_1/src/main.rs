use rayon::prelude::*;
use ndarray::{Array2, s};
// Adjust the imports to match your custom/forked rand version where rng() and random_range() exist
use rand::{rng, Rng};
use std::sync::{Arc, Mutex};
use std::fs::File;
use std::io::{Write, BufWriter};
use std::env;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut num_genes = 1000_usize;
    let mut num_samples = 50_usize;
    let mut output_file = String::from("partial_adjacency.bin");

    // Simple CLI argument parsing
    let args: Vec<String> = env::args().collect();
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--num-genes" => {
                i += 1;
                num_genes = args[i].parse::<usize>()?;
            }
            "--num-samples" => {
                i += 1;
                num_samples = args[i].parse::<usize>()?;
            }
            "--output" => {
                i += 1;
                output_file = args[i].clone();
            }
            _ => {}
        }
        i += 1;
    }

    println!("Number of genes: {}", num_genes);
    println!("Number of samples: {}", num_samples);
    println!("Output file: {}", output_file);

    // Synthetic data generation
    let expression_data = generate_synthetic_expression(num_genes, num_samples);

    // Shared adjacency matrix protected by a mutex
    let adjacency = Arc::new(Mutex::new(Array2::<f64>::zeros((num_genes, num_genes))));

    // Compute pairwise correlations in parallel (upper triangle only)
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

    // Obtain final matrix, write to binary file
    let final_adj_matrix = adjacency.lock().unwrap().clone();
    let file = File::create(&output_file)?;
    let mut writer = BufWriter::new(file);
    for i in 0..num_genes {
        for j in 0..num_genes {
            writer.write_all(&final_adj_matrix[[i, j]].to_ne_bytes())?;
        }
    }
    println!("Correlation adjacency matrix written to {}", output_file);
    Ok(())
}

// Generates synthetic gene expression data
fn generate_synthetic_expression(num_genes: usize, num_samples: usize) -> Array2<f64> {
    let mut data = Array2::<f64>::zeros((num_genes, num_samples));
    // Use the new rng() function
    let mut rng = rng(); 
    for i in 0..num_genes {
        for j in 0..num_samples {
            // Use the new random_range() function
            data[[i, j]] = rng.random_range(0.0..1000.0);
        }
    }
    data
}

// Computes Pearson correlation for two 1D slices
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
    if denom == 0.0 {
        0.0
    } else {
        num / denom
    }
}
