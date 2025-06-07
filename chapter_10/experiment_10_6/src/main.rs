use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use anyhow::{Context, Result};
use clap::Parser;
use log::{info, warn, debug};
use ndarray::{Array1, Array2};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use statrs::distribution::{Normal, Continuous};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input data file (CSV format)
    #[arg(short, long, default_value = "epigenomic_data.csv")]
    input: String,
    
    /// Output file for segmentation results
    #[arg(short, long, default_value = "segmentation_results.json")]
    output: String,
    
    /// Number of hidden states
    #[arg(short, long, default_value_t = 3)]
    states: usize,
    
    /// Maximum iterations for EM algorithm
    #[arg(short, long, default_value_t = 100)]
    max_iterations: usize,
    
    /// Convergence threshold
    #[arg(short, long, default_value_t = 1e-6)]
    tolerance: f64,
    
    /// Use parallel processing
    #[arg(long)]
    parallel: bool,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
struct StateParams {
    emission_mean: f64,
    emission_var: f64,
    transition_probs: Vec<f64>,
    initial_prob: f64,
}

#[derive(Serialize, Deserialize, Debug)]
struct HMMModel {
    states: Vec<StateParams>,
    num_states: usize,
    log_likelihood: f64,
    iterations: usize,
}

#[derive(Serialize, Deserialize, Debug)]
struct SegmentationResult {
    predicted_states: Vec<usize>,
    posterior_probs: Vec<Vec<f64>>,
    model: HMMModel,
    data_length: usize,
    processing_time_ms: u128,
}

#[derive(Debug)]
struct EpigenomicData {
    values: Vec<f64>,
    positions: Vec<usize>,
    chromosome: String,
}

fn main() -> Result<()> {
    env_logger::init();
    let args = Args::parse();
    
    info!("Starting epigenomic HMM segmentation");
    info!("Input file: {}", args.input);
    info!("Number of states: {}", args.states);
    info!("Parallel processing: {}", args.parallel);
    
    let start_time = std::time::Instant::now();
    
    // Load data
    let data = load_epigenomic_data(&args.input)
        .with_context(|| format!("Failed to load data from {}", args.input))?;
    
    info!("Loaded {} data points from chromosome {}", data.values.len(), data.chromosome);
    
    // Initialize HMM model
    let mut model = initialize_hmm_model(args.states, &data.values);
    
    // Train the model using EM algorithm
    let trained_model = if args.parallel {
        train_hmm_parallel(&data.values, model, args.max_iterations, args.tolerance)?
    } else {
        train_hmm(&data.values, model, args.max_iterations, args.tolerance)?
    };
    
    // Perform Viterbi decoding for state segmentation
    let (predicted_states, posterior_probs) = viterbi_decode(&data.values, &trained_model)?;
    
    let processing_time = start_time.elapsed().as_millis();
    
    // Create result structure
    let result = SegmentationResult {
        predicted_states,
        posterior_probs,
        model: trained_model,
        data_length: data.values.len(),
        processing_time_ms: processing_time,
    };
    
    // Save results
    save_results(&result, &args.output)
        .with_context(|| format!("Failed to save results to {}", args.output))?;
    
    info!("Segmentation completed in {} ms", processing_time);
    info!("Results saved to {}", args.output);
    
    // Print summary
    print_segmentation_summary(&result);
    
    Ok(())
}

fn load_epigenomic_data(filepath: &str) -> Result<EpigenomicData> {
    let file = File::open(filepath)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();
    
    // Skip header
    lines.next();
    
    let mut values = Vec::new();
    let mut positions = Vec::new();
    let mut chromosome = String::from("chr1");
    
    for (i, line) in lines.enumerate() {
        let line = line?;
        let parts: Vec<&str> = line.split(',').collect();
        
        if parts.len() >= 3 {
            if i == 0 {
                chromosome = parts[0].to_string();
            }
            let position: usize = parts[1].parse()
                .with_context(|| format!("Invalid position at line {}", i + 2))?;
            let value: f64 = parts[2].parse()
                .with_context(|| format!("Invalid value at line {}", i + 2))?;
            
            positions.push(position);
            values.push(value);
        }
    }
    
    if values.is_empty() {
        anyhow::bail!("No valid data found in input file");
    }
    
    Ok(EpigenomicData {
        values,
        positions,
        chromosome,
    })
}

fn initialize_hmm_model(num_states: usize, data: &[f64]) -> HMMModel {
    let mut rng = rand::thread_rng();
    
    // Calculate data statistics for initialization
    let data_mean = data.iter().sum::<f64>() / data.len() as f64;
    let data_var = data.iter()
        .map(|x| (x - data_mean).powi(2))
        .sum::<f64>() / data.len() as f64;
    
    let mut states = Vec::new();
    
    for i in 0..num_states {
        // Initialize emission parameters around data distribution
        let mean_offset = (i as f64 - (num_states as f64 - 1.0) / 2.0) * data_var.sqrt();
        let emission_mean = data_mean + mean_offset;
        let emission_var = data_var / num_states as f64;
        
        // Initialize transition probabilities (slightly favor self-transitions)
        let mut transition_probs = vec![0.1 / (num_states - 1) as f64; num_states];
        transition_probs[i] = 0.9;
        
        let initial_prob = 1.0 / num_states as f64;
        
        states.push(StateParams {
            emission_mean,
            emission_var,
            transition_probs,
            initial_prob,
        });
    }
    
    HMMModel {
        states,
        num_states,
        log_likelihood: f64::NEG_INFINITY,
        iterations: 0,
    }
}

fn train_hmm(data: &[f64], mut model: HMMModel, max_iterations: usize, tolerance: f64) -> Result<HMMModel> {
    let n = data.len();
    let k = model.num_states;
    
    for iteration in 0..max_iterations {
        debug!("EM iteration {}", iteration + 1);
        
        // E-step: Forward-backward algorithm
        let (alpha, beta, log_likelihood) = forward_backward(data, &model)?;
        
        // Check convergence
        let improvement = log_likelihood - model.log_likelihood;
        if iteration > 0 && improvement.abs() < tolerance {
            info!("Converged after {} iterations", iteration + 1);
            model.log_likelihood = log_likelihood;
            model.iterations = iteration + 1;
            break;
        }
        
        model.log_likelihood = log_likelihood;
        
        // M-step: Update parameters
        let gamma = compute_gamma(&alpha, &beta)?;
        let xi = compute_xi(data, &model, &alpha, &beta)?;
        
        // Update initial probabilities
        for i in 0..k {
            model.states[i].initial_prob = gamma[[0, i]];
        }
        
        // Update transition probabilities
        for i in 0..k {
            let xi_sum: f64 = (0..n-1).map(|t| (0..k).map(|j| xi[[t, i, j]]).sum::<f64>()).sum();
            for j in 0..k {
                let xi_ij_sum: f64 = (0..n-1).map(|t| xi[[t, i, j]]).sum();
                model.states[i].transition_probs[j] = if xi_sum > 0.0 { xi_ij_sum / xi_sum } else { 1.0 / k as f64 };
            }
        }
        
        // Update emission parameters
        for i in 0..k {
            let gamma_sum: f64 = (0..n).map(|t| gamma[[t, i]]).sum();
            
            if gamma_sum > 0.0 {
                // Update mean
                let weighted_sum: f64 = (0..n).map(|t| gamma[[t, i]] * data[t]).sum();
                model.states[i].emission_mean = weighted_sum / gamma_sum;
                
                // Update variance
                let variance_sum: f64 = (0..n).map(|t| {
                    gamma[[t, i]] * (data[t] - model.states[i].emission_mean).powi(2)
                }).sum();
                model.states[i].emission_var = (variance_sum / gamma_sum).max(1e-6);
            }
        }
        
        info!("Iteration {}: Log-likelihood = {:.6}", iteration + 1, log_likelihood);
    }
    
    model.iterations = max_iterations.min(model.iterations);
    Ok(model)
}

fn train_hmm_parallel(data: &[f64], model: HMMModel, max_iterations: usize, tolerance: f64) -> Result<HMMModel> {
    // For simplicity, this uses the same algorithm but could be optimized for parallel processing
    info!("Using parallel processing (simplified implementation)");
    train_hmm(data, model, max_iterations, tolerance)
}

fn forward_backward(data: &[f64], model: &HMMModel) -> Result<(Array2<f64>, Array2<f64>, f64)> {
    let n = data.len();
    let k = model.num_states;
    
    let mut alpha = Array2::zeros((n, k));
    let mut beta = Array2::zeros((n, k));
    
    // Forward pass
    for i in 0..k {
        let emission_prob = emission_probability(data[0], &model.states[i]);
        alpha[[0, i]] = model.states[i].initial_prob.ln() + emission_prob.ln();
    }
    
    for t in 1..n {
        for j in 0..k {
            let emission_prob = emission_probability(data[t], &model.states[j]);
            let mut sum = f64::NEG_INFINITY;
            
            for i in 0..k {
                let trans_prob = model.states[i].transition_probs[j];
                if trans_prob > 0.0 {
                    let prob = alpha[[t-1, i]] + trans_prob.ln();
                    sum = log_sum_exp(sum, prob);
                }
            }
            
            alpha[[t, j]] = sum + emission_prob.ln();
        }
    }
    
    // Calculate log-likelihood
    let log_likelihood = (0..k)
        .map(|i| alpha[[n-1, i]])
        .fold(f64::NEG_INFINITY, log_sum_exp);
    
    // Backward pass
    for i in 0..k {
        beta[[n-1, i]] = 0.0; // log(1) = 0
    }
    
    for t in (0..n-1).rev() {
        for i in 0..k {
            let mut sum = f64::NEG_INFINITY;
            
            for j in 0..k {
                let trans_prob = model.states[i].transition_probs[j];
                if trans_prob > 0.0 {
                    let emission_prob = emission_probability(data[t+1], &model.states[j]);
                    let prob = trans_prob.ln() + emission_prob.ln() + beta[[t+1, j]];
                    sum = log_sum_exp(sum, prob);
                }
            }
            
            beta[[t, i]] = sum;
        }
    }
    
    Ok((alpha, beta, log_likelihood))
}

fn compute_gamma(alpha: &Array2<f64>, beta: &Array2<f64>) -> Result<Array2<f64>> {
    let (n, k) = alpha.dim();
    let mut gamma = Array2::zeros((n, k));
    
    for t in 0..n {
        let log_norm = (0..k)
            .map(|i| alpha[[t, i]] + beta[[t, i]])
            .fold(f64::NEG_INFINITY, log_sum_exp);
        
        for i in 0..k {
            gamma[[t, i]] = (alpha[[t, i]] + beta[[t, i]] - log_norm).exp();
        }
    }
    
    Ok(gamma)
}

fn compute_xi(data: &[f64], model: &HMMModel, alpha: &Array2<f64>, beta: &Array2<f64>) -> Result<ndarray::Array3<f64>> {
    let n = data.len();
    let k = model.num_states;
    let mut xi = ndarray::Array3::zeros((n-1, k, k));
    
    for t in 0..n-1 {
        let mut log_norm = f64::NEG_INFINITY;
        
        // Calculate normalization constant
        for i in 0..k {
            for j in 0..k {
                let trans_prob = model.states[i].transition_probs[j];
                if trans_prob > 0.0 {
                    let emission_prob = emission_probability(data[t+1], &model.states[j]);
                    let prob = alpha[[t, i]] + trans_prob.ln() + emission_prob.ln() + beta[[t+1, j]];
                    log_norm = log_sum_exp(log_norm, prob);
                }
            }
        }
        
        // Calculate xi
        for i in 0..k {
            for j in 0..k {
                let trans_prob = model.states[i].transition_probs[j];
                if trans_prob > 0.0 {
                    let emission_prob = emission_probability(data[t+1], &model.states[j]);
                    let log_prob = alpha[[t, i]] + trans_prob.ln() + emission_prob.ln() + beta[[t+1, j]] - log_norm;
                    xi[[t, i, j]] = log_prob.exp();
                } else {
                    xi[[t, i, j]] = 0.0;
                }
            }
        }
    }
    
    Ok(xi)
}

fn viterbi_decode(data: &[f64], model: &HMMModel) -> Result<(Vec<usize>, Vec<Vec<f64>>)> {
    let n = data.len();
    let k = model.num_states;
    
    let mut delta = Array2::zeros((n, k));
    let mut psi = Array2::zeros((n, k));
    
    // Initialize
    for i in 0..k {
        let emission_prob = emission_probability(data[0], &model.states[i]);
        delta[[0, i]] = model.states[i].initial_prob.ln() + emission_prob.ln();
        psi[[0, i]] = 0.0;
    }
    
    // Forward pass
    for t in 1..n {
        for j in 0..k {
            let emission_prob = emission_probability(data[t], &model.states[j]);
            let mut max_prob = f64::NEG_INFINITY;
            let mut best_prev = 0;
            
            for i in 0..k {
                let trans_prob = model.states[i].transition_probs[j];
                if trans_prob > 0.0 {
                    let prob = delta[[t-1, i]] + trans_prob.ln();
                    if prob > max_prob {
                        max_prob = prob;
                        best_prev = i;
                    }
                }
            }
            
            delta[[t, j]] = max_prob + emission_prob.ln();
            psi[[t, j]] = best_prev as f64;
        }
    }
    
    // Find best final state
    let mut best_final_state = 0;
    let mut max_final_prob = f64::NEG_INFINITY;
    for i in 0..k {
        if delta[[n-1, i]] > max_final_prob {
            max_final_prob = delta[[n-1, i]];
            best_final_state = i;
        }
    }
    
    // Backward pass (traceback)
    let mut path = vec![0; n];
    path[n-1] = best_final_state;
    
    for t in (1..n).rev() {
        path[t-1] = psi[[t, path[t]]] as usize;
    }
    
    // Compute posterior probabilities using forward-backward
    let (alpha, beta, _) = forward_backward(data, model)?;
    let gamma = compute_gamma(&alpha, &beta)?;
    
    let mut posterior_probs = Vec::new();
    for t in 0..n {
        let mut probs = Vec::new();
        for i in 0..k {
            probs.push(gamma[[t, i]]);
        }
        posterior_probs.push(probs);
    }
    
    Ok((path, posterior_probs))
}

fn emission_probability(observation: f64, state: &StateParams) -> f64 {
    let normal = Normal::new(state.emission_mean, state.emission_var.sqrt()).unwrap();
    normal.pdf(observation).max(1e-300) // Avoid underflow
}

fn log_sum_exp(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY && b == f64::NEG_INFINITY {
        f64::NEG_INFINITY
    } else if a > b {
        a + (b - a).exp().ln_1p()
    } else {
        b + (a - b).exp().ln_1p()
    }
}

fn save_results(result: &SegmentationResult, filepath: &str) -> Result<()> {
    let json = serde_json::to_string_pretty(result)?;
    let mut file = File::create(filepath)?;
    file.write_all(json.as_bytes())?;
    Ok(())
}

fn print_segmentation_summary(result: &SegmentationResult) {
    println!("\n=== Segmentation Summary ===");
    println!("Data length: {}", result.data_length);
    println!("Number of states: {}", result.model.num_states);
    println!("Final log-likelihood: {:.6}", result.model.log_likelihood);
    println!("Training iterations: {}", result.model.iterations);
    println!("Processing time: {} ms", result.processing_time_ms);
    
    // State distribution
    let mut state_counts = vec![0; result.model.num_states];
    for &state in &result.predicted_states {
        state_counts[state] += 1;
    }
    
    println!("\nState distribution:");
    for (i, count) in state_counts.iter().enumerate() {
        let percentage = *count as f64 / result.data_length as f64 * 100.0;
        println!("  State {}: {} positions ({:.1}%)", i, count, percentage);
    }
    
    println!("\nModel parameters:");
    for (i, state) in result.model.states.iter().enumerate() {
        println!("  State {}: mean={:.3}, var={:.3}, initial_prob={:.3}", 
                 i, state.emission_mean, state.emission_var, state.initial_prob);
    }
}