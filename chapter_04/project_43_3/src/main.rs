use rand::{rng, Rng};
use rayon; // We'll call `rayon::scope()` directly.
use std::sync::{Arc, Mutex};

/// A structure to hold the state of a Gibbs sampling chain:
/// - motif_positions: the current motif start positions in each sequence
/// - k: motif length
/// - sequences: shared references to the input data
#[derive(Debug)] // <-- Add Debug to fix the E0277 error
pub struct GibbsSampler {
    motif_positions: Vec<usize>,
    k: usize,
    sequences: Vec<Vec<u8>>,
}

impl GibbsSampler {
    /// Creates a new chain with randomly initialized motif positions.
    pub fn new(sequences: Vec<Vec<u8>>, k: usize) -> Self {
        // Replace `thread_rng()` with `rng()`
        let mut rng = rng();
        let mut motif_positions = Vec::new();

        for seq in &sequences {
            // random starting position, ensuring space for the motif
            if seq.len() >= k {
                // Replace `gen_range(0..=x)` with `random_range(0..=x)`
                let start = rng.random_range(0..=(seq.len() - k));
                motif_positions.push(start);
            } else {
                motif_positions.push(0); // invalid, but for demonstration
            }
        }

        GibbsSampler {
            motif_positions,
            k,
            sequences,
        }
    }

    /// Runs one iteration of Gibbs sampling:
    /// picks a random sequence, "unassigns" its motif, then samples a new position from a distribution.
    pub fn run_one_iteration(&mut self) {
        let mut rng = rng(); // again, replace if your environment requires
        let seq_index = rng.random_range(0..self.sequences.len());
        let seq = &self.sequences[seq_index];

        // Remove seq_index from motif model (in real code, recalc partial counts or re-estimate PWM)
        // Then sample new position based on a simplistic probability model:
        let mut cumulative_probs = Vec::new();
        let mut sum_probs = 0.0;

        for start in 0..=(seq.len().saturating_sub(self.k)) {
            // compute some pseudo-prob (like counting 'A', or any placeholder measure)
            let prob = self.motif_likelihood(seq, start);
            sum_probs += prob;
            cumulative_probs.push(sum_probs);
        }

        let r = rng.random_range(0.0..sum_probs);
        // find where r fits into the cumulative distribution
        let pos = match cumulative_probs.iter().position(|&cp| cp >= r) {
            Some(idx) => idx,
            None => 0,
        };
        self.motif_positions[seq_index] = pos;
    }

    /// A toy likelihood function for demonstration: simply counts A's in the k-mer
    /// more A's => higher probability
    fn motif_likelihood(&self, seq: &[u8], start: usize) -> f64 {
        let mut score = 1.0;
        let end = start + self.k;
        for i in start..end {
            if i < seq.len() && (seq[i] == b'A' || seq[i] == b'a') {
                score += 1.0;
            }
        }
        score
    }
}

/// Runs multiple chains of Gibbs sampling in parallel using rayon's concurrency.
pub fn run_parallel_chains(
    sequences: Vec<Vec<u8>>,
    k: usize,
    num_chains: usize,
    iterations: usize,
) -> Vec<GibbsSampler> {
    let sampler_results = Arc::new(Mutex::new(Vec::new()));

    rayon::scope(|s| {
        for _ in 0..num_chains {
            let s_res = Arc::clone(&sampler_results);
            let seq_clone = sequences.clone();
            s.spawn(move |_| {
                let mut sampler = GibbsSampler::new(seq_clone, k);
                for _iter in 0..iterations {
                    sampler.run_one_iteration();
                }
                // store final sampler in results
                let mut lock = s_res.lock().unwrap();
                lock.push(sampler);
            });
        }
    });

    // Unwrap the Arc/Mutex to get the final results
    let final_samplers = Arc::try_unwrap(sampler_results)
        .unwrap()
        .into_inner()
        .unwrap();
    final_samplers
}

fn main() {
    // Suppose we load sequences from a file or HPC ephemeral container; here, just a few toy reads:
    let sequences = vec![
        b"ACGATGATGAC".to_vec(),
        b"TTTTAAAACCCCGG".to_vec(),
        b"AAAATGATGAAAA".to_vec(),
    ];
    let k = 5;
    let num_chains = 3;
    let iterations = 10;

    // run multiple chains in parallel
    let samplers = run_parallel_chains(sequences, k, num_chains, iterations);

    // In real HPC usage, we might combine results from these samplers (like picking the best chain).
    for (i, sampler) in samplers.iter().enumerate() {
        println!("Sampler {} final motif positions: {:?}", i, sampler.motif_positions);
    }
}
