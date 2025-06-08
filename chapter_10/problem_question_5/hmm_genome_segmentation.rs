use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ChromatinState {
    Active = 0,
    Repressed = 1,
    Poised = 2,
}

impl ChromatinState {
    pub fn from_index(index: usize) -> Self {
        match index {
            0 => ChromatinState::Active,
            1 => ChromatinState::Repressed,
            2 => ChromatinState::Poised,
            _ => panic!("Invalid chromatin state index"),
        }
    }

    pub fn to_string(&self) -> &'static str {
        match self {
            ChromatinState::Active => "Active",
            ChromatinState::Repressed => "Repressed",
            ChromatinState::Poised => "Poised",
        }
    }
}

#[derive(Debug, Clone)]
pub struct HistoneData {
    pub h3k4me3: f64,  // Active promoter mark
    pub h3k27ac: f64,  // Active enhancer mark
    pub h3k27me3: f64, // Repressive mark (polycomb)
    pub h3k4me1: f64,  // Enhancer mark
}

impl HistoneData {
    pub fn new(h3k4me3: f64, h3k27ac: f64, h3k27me3: f64, h3k4me1: f64) -> Self {
        Self { h3k4me3, h3k27ac, h3k27me3, h3k4me1 }
    }

    pub fn to_vector(&self) -> Vec<f64> {
        vec![self.h3k4me3, self.h3k27ac, self.h3k27me3, self.h3k4me1]
    }
}

#[derive(Debug, Clone)]
pub struct HMMParameters {
    pub num_states: usize,
    pub num_features: usize,
    pub transition_matrix: Vec<Vec<f64>>,
    pub emission_means: Vec<Vec<f64>>,
    pub emission_stds: Vec<Vec<f64>>,
    pub initial_probs: Vec<f64>,
}

impl HMMParameters {
    pub fn new(num_states: usize, num_features: usize) -> Self {
        // Initialize with reasonable biological priors
        let mut transition_matrix = vec![vec![0.0; num_states]; num_states];
        let mut emission_means = vec![vec![0.0; num_features]; num_states];
        let mut emission_stds = vec![vec![1.0; num_features]; num_states];
        let initial_probs = vec![1.0 / num_states as f64; num_states];

        // Set biologically informed initial parameters
        // State 0: Active (high H3K4me3, H3K27ac, low H3K27me3)
        emission_means[0] = vec![3.0, 2.5, 0.5, 1.5]; // H3K4me3, H3K27ac, H3K27me3, H3K4me1
        
        // State 1: Repressed (low active marks, high H3K27me3)
        emission_means[1] = vec![0.3, 0.2, 3.0, 0.4];
        
        // State 2: Poised (moderate H3K4me3, high H3K27me3 - bivalent)
        emission_means[2] = vec![2.0, 0.5, 2.5, 1.0];

        // Initialize transition probabilities with self-transition bias
        for i in 0..num_states {
            for j in 0..num_states {
                if i == j {
                    transition_matrix[i][j] = 0.8; // Strong self-transition
                } else {
                    transition_matrix[i][j] = 0.1; // Weak transitions to other states
                }
            }
        }

        Self {
            num_states,
            num_features,
            transition_matrix,
            emission_means,
            emission_stds,
            initial_probs,
        }
    }

    pub fn normalize_transitions(&mut self) {
        for i in 0..self.num_states {
            let sum: f64 = self.transition_matrix[i].iter().sum();
            if sum > 0.0 {
                for j in 0..self.num_states {
                    self.transition_matrix[i][j] /= sum;
                }
            }
        }
    }
}

pub struct GenomeHMM {
    pub params: HMMParameters,
    pub convergence_threshold: f64,
    pub max_iterations: usize,
}

impl GenomeHMM {
    pub fn new(num_states: usize, num_features: usize) -> Self {
        Self {
            params: HMMParameters::new(num_states, num_features),
            convergence_threshold: 1e-6,
            max_iterations: 100,
        }
    }

    // Gaussian emission probability
    fn gaussian_prob(&self, value: f64, mean: f64, std: f64) -> f64 {
        let exp_part = -0.5 * ((value - mean) / std).powi(2);
        (1.0 / (std * (2.0 * std::f64::consts::PI).sqrt())) * exp_part.exp()
    }

    // Calculate emission probability for multivariate observation
    fn emission_prob(&self, state: usize, observation: &[f64]) -> f64 {
        let mut prob = 1.0;
        for (i, &value) in observation.iter().enumerate() {
            let mean = self.params.emission_means[state][i];
            let std = self.params.emission_stds[state][i];
            prob *= self.gaussian_prob(value, mean, std);
        }
        prob.max(1e-300) // Prevent underflow
    }

    // Viterbi algorithm for finding most likely state sequence
    pub fn viterbi(&self, observations: &[HistoneData]) -> Vec<ChromatinState> {
        let n_obs = observations.len();
        let n_states = self.params.num_states;
        
        if n_obs == 0 {
            return Vec::new();
        }

        // Initialize Viterbi tables
        let mut viterbi_table = vec![vec![0.0; n_states]; n_obs];
        let mut path_table = vec![vec![0; n_states]; n_obs];

        // Initialize first observation
        for state in 0..n_states {
            let obs_vec = observations[0].to_vector();
            viterbi_table[0][state] = 
                self.params.initial_probs[state].ln() + 
                self.emission_prob(state, &obs_vec).ln();
        }

        // Forward pass
        for t in 1..n_obs {
            let obs_vec = observations[t].to_vector();
            for curr_state in 0..n_states {
                let mut best_prob = f64::NEG_INFINITY;
                let mut best_prev_state = 0;

                for prev_state in 0..n_states {
                    let prob = viterbi_table[t - 1][prev_state] +
                              self.params.transition_matrix[prev_state][curr_state].ln() +
                              self.emission_prob(curr_state, &obs_vec).ln();
                    
                    if prob > best_prob {
                        best_prob = prob;
                        best_prev_state = prev_state;
                    }
                }

                viterbi_table[t][curr_state] = best_prob;
                path_table[t][curr_state] = best_prev_state;
            }
        }

        // Find best final state
        let mut best_final_state = 0;
        let mut best_final_prob = viterbi_table[n_obs - 1][0];
        for state in 1..n_states {
            if viterbi_table[n_obs - 1][state] > best_final_prob {
                best_final_prob = viterbi_table[n_obs - 1][state];
                best_final_state = state;
            }
        }

        // Backtrack to find optimal path
        let mut path = vec![0; n_obs];
        path[n_obs - 1] = best_final_state;
        for t in (1..n_obs).rev() {
            path[t - 1] = path_table[t][path[t]];
        }

        path.into_iter().map(ChromatinState::from_index).collect()
    }

    // Forward-backward algorithm for parameter estimation
    pub fn forward_backward(&self, observations: &[HistoneData]) -> (Vec<Vec<f64>>, f64) {
        let n_obs = observations.len();
        let n_states = self.params.num_states;
        
        if n_obs == 0 {
            return (Vec::new(), 0.0);
        }

        // Forward pass
        let mut forward = vec![vec![0.0; n_states]; n_obs];
        let mut scaling_factors = vec![0.0; n_obs];

        // Initialize forward probabilities
        let obs_vec = observations[0].to_vector();
        for state in 0..n_states {
            forward[0][state] = self.params.initial_probs[state] * 
                               self.emission_prob(state, &obs_vec);
        }
        
        // Scale to prevent underflow
        let sum: f64 = forward[0].iter().sum();
        scaling_factors[0] = sum;
        if sum > 0.0 {
            for state in 0..n_states {
                forward[0][state] /= sum;
            }
        }

        // Forward recursion
        for t in 1..n_obs {
            let obs_vec = observations[t].to_vector();
            for curr_state in 0..n_states {
                forward[t][curr_state] = 0.0;
                for prev_state in 0..n_states {
                    forward[t][curr_state] += forward[t - 1][prev_state] *
                                            self.params.transition_matrix[prev_state][curr_state];
                }
                forward[t][curr_state] *= self.emission_prob(curr_state, &obs_vec);
            }
            
            // Scale
            let sum: f64 = forward[t].iter().sum();
            scaling_factors[t] = sum;
            if sum > 0.0 {
                for state in 0..n_states {
                    forward[t][state] /= sum;
                }
            }
        }

        // Backward pass
        let mut backward = vec![vec![0.0; n_states]; n_obs];
        for state in 0..n_states {
            backward[n_obs - 1][state] = 1.0;
        }

        for t in (0..n_obs - 1).rev() {
            let next_obs_vec = observations[t + 1].to_vector();
            for curr_state in 0..n_states {
                backward[t][curr_state] = 0.0;
                for next_state in 0..n_states {
                    backward[t][curr_state] += 
                        self.params.transition_matrix[curr_state][next_state] *
                        self.emission_prob(next_state, &next_obs_vec) *
                        backward[t + 1][next_state];
                }
                backward[t][curr_state] /= scaling_factors[t + 1];
            }
        }

        // Calculate gamma (posterior probabilities)
        let mut gamma = vec![vec![0.0; n_states]; n_obs];
        for t in 0..n_obs {
            let sum: f64 = (0..n_states)
                .map(|s| forward[t][s] * backward[t][s])
                .sum();
            
            if sum > 0.0 {
                for state in 0..n_states {
                    gamma[t][state] = (forward[t][state] * backward[t][state]) / sum;
                }
            }
        }

        // Calculate log likelihood
        let log_likelihood: f64 = scaling_factors.iter()
            .map(|&s| if s > 0.0 { s.ln() } else { 0.0 })
            .sum();

        (gamma, log_likelihood)
    }

    // Baum-Welch algorithm for training
    pub fn train(&mut self, observations: &[HistoneData]) -> f64 {
        let mut prev_log_likelihood = f64::NEG_INFINITY;
        
        for iteration in 0..self.max_iterations {
            let (gamma, log_likelihood) = self.forward_backward(observations);
            
            println!("Iteration {}: Log-likelihood = {:.6}", iteration + 1, log_likelihood);
            
            // Check convergence
            if (log_likelihood - prev_log_likelihood).abs() < self.convergence_threshold {
                println!("Converged after {} iterations", iteration + 1);
                break;
            }
            
            prev_log_likelihood = log_likelihood;
            
            // Update parameters
            self.update_parameters(observations, &gamma);
        }
        
        prev_log_likelihood
    }

    fn update_parameters(&mut self, observations: &[HistoneData], gamma: &[Vec<f64>]) {
        let n_obs = observations.len();
        let n_states = self.params.num_states;
        let n_features = self.params.num_features;

        // Update initial probabilities
        for state in 0..n_states {
            self.params.initial_probs[state] = gamma[0][state];
        }

        // Update emission parameters
        for state in 0..n_states {
            let gamma_sum: f64 = gamma.iter().map(|g| g[state]).sum();
            
            if gamma_sum > 0.0 {
                // Update means
                for feature in 0..n_features {
                    let weighted_sum: f64 = observations.iter()
                        .zip(gamma.iter())
                        .map(|(obs, g)| obs.to_vector()[feature] * g[state])
                        .sum();
                    self.params.emission_means[state][feature] = weighted_sum / gamma_sum;
                }

                // Update standard deviations
                for feature in 0..n_features {
                    let mean = self.params.emission_means[state][feature];
                    let weighted_var_sum: f64 = observations.iter()
                        .zip(gamma.iter())
                        .map(|(obs, g)| {
                            let diff = obs.to_vector()[feature] - mean;
                            diff * diff * g[state]
                        })
                        .sum();
                    self.params.emission_stds[state][feature] = 
                        (weighted_var_sum / gamma_sum).sqrt().max(0.1);
                }
            }
        }

        // Update transition probabilities
        for i in 0..n_states {
            let mut transition_sum = vec![0.0; n_states];
            let mut total_transitions = 0.0;

            for t in 0..n_obs - 1 {
                let next_obs = observations[t + 1].to_vector();
                let normalizer: f64 = (0..n_states)
                    .map(|j| gamma[t][i] * 
                           self.params.transition_matrix[i][j] * 
                           self.emission_prob(j, &next_obs))
                    .sum();

                if normalizer > 0.0 {
                    for j in 0..n_states {
                        let xi = (gamma[t][i] * 
                                 self.params.transition_matrix[i][j] * 
                                 self.emission_prob(j, &next_obs)) / normalizer;
                        transition_sum[j] += xi;
                        total_transitions += xi;
                    }
                }
            }

            if total_transitions > 0.0 {
                for j in 0..n_states {
                    self.params.transition_matrix[i][j] = transition_sum[j] / total_transitions;
                }
            }
        }

        self.params.normalize_transitions();
    }

    pub fn print_parameters(&self) {
        println!("\n=== HMM Parameters ===");
        
        println!("\nInitial Probabilities:");
        for (i, &prob) in self.params.initial_probs.iter().enumerate() {
            println!("  {}: {:.4}", ChromatinState::from_index(i).to_string(), prob);
        }
        
        println!("\nTransition Matrix:");
        print!("        ");
        for j in 0..self.params.num_states {
            print!("{:>12}", ChromatinState::from_index(j).to_string());
        }
        println!();
        
        for i in 0..self.params.num_states {
            print!("{:>8}", ChromatinState::from_index(i).to_string());
            for j in 0..self.params.num_states {
                print!("{:>12.4}", self.params.transition_matrix[i][j]);
            }
            println!();
        }
        
        println!("\nEmission Parameters:");
        let feature_names = ["H3K4me3", "H3K27ac", "H3K27me3", "H3K4me1"];
        for (state_idx, state_name) in [
            ChromatinState::Active.to_string(),
            ChromatinState::Repressed.to_string(),
            ChromatinState::Poised.to_string()
        ].iter().enumerate() {
            println!("  {}:", state_name);
            for (feature_idx, &feature_name) in feature_names.iter().enumerate() {
                println!("    {}: μ={:.3}, σ={:.3}", 
                        feature_name,
                        self.params.emission_means[state_idx][feature_idx],
                        self.params.emission_stds[state_idx][feature_idx]);
            }
        }
    }
}

// Data generation and file I/O functions
pub fn generate_synthetic_data(length: usize, seed: u64) -> Vec<HistoneData> {
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};
    
    let mut rng_state = seed;
    let mut next_random = || {
        rng_state = rng_state.wrapping_mul(1103515245).wrapping_add(12345);
        (rng_state >> 16) as f64 / 32768.0
    };

    let mut data = Vec::with_capacity(length);
    let mut current_state = 0; // Start in Active state
    
    for _ in 0..length {
        // Occasional state transitions (simulate domain structure)
        if next_random() < 0.05 {
            current_state = (next_random() * 3.0) as usize;
        }
        
        let (h3k4me3, h3k27ac, h3k27me3, h3k4me1) = match current_state {
            0 => { // Active
                let noise_scale = 0.5;
                (3.0 + (next_random() - 0.5) * noise_scale,
                 2.5 + (next_random() - 0.5) * noise_scale,
                 0.5 + (next_random() - 0.5) * noise_scale,
                 1.5 + (next_random() - 0.5) * noise_scale)
            },
            1 => { // Repressed
                let noise_scale = 0.3;
                (0.3 + (next_random() - 0.5) * noise_scale,
                 0.2 + (next_random() - 0.5) * noise_scale,
                 3.0 + (next_random() - 0.5) * noise_scale,
                 0.4 + (next_random() - 0.5) * noise_scale)
            },
            _ => { // Poised
                let noise_scale = 0.4;
                (2.0 + (next_random() - 0.5) * noise_scale,
                 0.5 + (next_random() - 0.5) * noise_scale,
                 2.5 + (next_random() - 0.5) * noise_scale,
                 1.0 + (next_random() - 0.5) * noise_scale)
            }
        };
        
        data.push(HistoneData::new(
            h3k4me3.max(0.0),
            h3k27ac.max(0.0), 
            h3k27me3.max(0.0),
            h3k4me1.max(0.0)
        ));
    }
    
    data
}

pub fn write_segmentation_results(
    states: &[ChromatinState], 
    filename: &str
) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = File::create(filename)?;
    writeln!(file, "position\tstate\tstate_name")?;
    
    for (i, &state) in states.iter().enumerate() {
        writeln!(file, "{}\t{}\t{}", i + 1, state as usize, state.to_string())?;
    }
    
    Ok(())
}

pub fn evaluate_segmentation_accuracy(
    predicted: &[ChromatinState],
    true_states: &[ChromatinState]
) -> f64 {
    if predicted.len() != true_states.len() {
        return 0.0;
    }
    
    let correct = predicted.iter()
        .zip(true_states.iter())
        .filter(|(&p, &t)| p == t)
        .count();
    
    correct as f64 / predicted.len() as f64
}

// Main function demonstrating usage
fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== HMM Genome Segmentation Demo ===\n");
    
    // Generate synthetic data
    println!("Generating synthetic genomic data...");
    let training_data = generate_synthetic_data(1000, 42);
    let test_data = generate_synthetic_data(500, 123);
    
    // Initialize and train HMM
    println!("Initializing HMM with 3 states and 4 histone marks...");
    let mut hmm = GenomeHMM::new(3, 4);
    
    println!("\nInitial parameters:");
    hmm.print_parameters();
    
    println!("\nTraining HMM...");
    let final_likelihood = hmm.train(&training_data);
    
    println!("\nFinal parameters after training:");
    hmm.print_parameters();
    
    // Perform segmentation using Viterbi algorithm
    println!("\nPerforming genome segmentation...");
    let predicted_states = hmm.viterbi(&test_data);
    
    // Generate true states for evaluation (simplified)
    let true_states = generate_true_states_for_evaluation(&test_data);
    
    // Evaluate accuracy
    let accuracy = evaluate_segmentation_accuracy(&predicted_states, &true_states);
    println!("Segmentation accuracy: {:.2}%", accuracy * 100.0);
    
    // Write results to file
    write_segmentation_results(&predicted_states, "genome_segmentation.txt")?;
    println!("Segmentation results written to 'genome_segmentation.txt'");
    
    // Print state distribution
    let mut state_counts = [0; 3];
    for &state in &predicted_states {
        state_counts[state as usize] += 1;
    }
    
    println!("\nState distribution in segmented genome:");
    for (i, &count) in state_counts.iter().enumerate() {
        let percentage = count as f64 / predicted_states.len() as f64 * 100.0;
        println!("  {}: {} positions ({:.1}%)", 
                ChromatinState::from_index(i).to_string(), count, percentage);
    }
    
    Ok(())
}

// Helper function to generate true states for evaluation
fn generate_true_states_for_evaluation(data: &[HistoneData]) -> Vec<ChromatinState> {
    data.iter().map(|obs| {
        // Simple rule-based classification for evaluation
        if obs.h3k4me3 > 2.0 && obs.h3k27ac > 1.5 && obs.h3k27me3 < 1.5 {
            ChromatinState::Active
        } else if obs.h3k27me3 > 2.0 && obs.h3k4me3 < 1.0 {
            ChromatinState::Repressed
        } else {
            ChromatinState::Poised
        }
    }).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hmm_initialization() {
        let hmm = GenomeHMM::new(3, 4);
        assert_eq!(hmm.params.num_states, 3);
        assert_eq!(hmm.params.num_features, 4);
        assert_eq!(hmm.params.transition_matrix.len(), 3);
        assert_eq!(hmm.params.emission_means.len(), 3);
    }
    
    #[test]
    fn test_viterbi_basic() {
        let hmm = GenomeHMM::new(3, 4);
        let data = vec![
            HistoneData::new(3.0, 2.5, 0.5, 1.5), // Should be Active
            HistoneData::new(0.3, 0.2, 3.0, 0.4), // Should be Repressed
        ];
        
        let states = hmm.viterbi(&data);
        assert_eq!(states.len(), 2);
    }
    
    #[test]
    fn test_synthetic_data_generation() {
        let data = generate_synthetic_data(100, 42);
        assert_eq!(data.len(), 100);
        
        // Check that all values are non-negative
        for obs in &data {
            assert!(obs.h3k4me3 >= 0.0);
            assert!(obs.h3k27ac >= 0.0);
            assert!(obs.h3k27me3 >= 0.0);
            assert!(obs.h3k4me1 >= 0.0);
        }
    }
    
    #[test]
    fn test_accuracy_calculation() {
        let predicted = vec![ChromatinState::Active, ChromatinState::Repressed];
        let truth = vec![ChromatinState::Active, ChromatinState::Active];
        let accuracy = evaluate_segmentation_accuracy(&predicted, &truth);
        assert_eq!(accuracy, 0.5);
    }
}