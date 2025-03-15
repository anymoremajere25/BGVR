use bio::io::fastq;
use std::env;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::collections::HashMap;

/// A simple function to build a Position Weight Matrix from FASTQ sequences.
/// We assume all sequences have the same length for simplicity.
fn build_pwm(sequences: &Vec<String>) -> Vec<HashMap<char, f64>> {
    // We expect sequences to be the same length. For real data, handle variable lengths carefully.
    let seq_length = sequences[0].len();
    let mut counts_per_position = vec![HashMap::new(); seq_length];

    // Initialize counts
    for position in 0..seq_length {
        for base in &['A', 'C', 'G', 'T'] {
            counts_per_position[position].insert(*base, 0.0);
        }
    }

    // Count nucleotides per position
    for seq in sequences {
        for (i, base) in seq.chars().enumerate() {
            if let Some(count) = counts_per_position[i].get_mut(&base) {
                *count += 1.0;
            }
        }
    }

    // Convert counts to probabilities (position-wise normalization)
    for position in 0..seq_length {
        let total_count: f64 = counts_per_position[position].values().sum();
        for base in &['A', 'C', 'G', 'T'] {
            let count = counts_per_position[position][base];
            let probability = if total_count > 0.0 {
                count / total_count
            } else {
                0.0
            };
            counts_per_position[position].insert(*base, probability);
        }
    }

    counts_per_position
}

/// A simple function to construct a Markov Random Field from the same sequences.
/// For demonstration, we use a 1st-order Markov chain, counting transitions between consecutive bases.
fn build_mrf(sequences: &Vec<String>) -> HashMap<(char, char), f64> {
    // Transition counts for pairs (X->Y)
    let mut transition_counts = HashMap::new();
    for &base1 in &['A', 'C', 'G', 'T'] {
        for &base2 in &['A', 'C', 'G', 'T'] {
            transition_counts.insert((base1, base2), 0.0);
        }
    }

    let mut total_pairs = 0.0;
    for seq in sequences {
        let chars: Vec<char> = seq.chars().collect();
        for i in 0..chars.len()-1 {
            if let Some(count) = transition_counts.get_mut(&(chars[i], chars[i+1])) {
                *count += 1.0;
                total_pairs += 1.0;
            }
        }
    }

    // Convert counts to transition probabilities
    let mut transition_probabilities = HashMap::new();
    for (&pair, &count) in &transition_counts {
        let prob = if total_pairs > 0.0 {
            count / total_pairs
        } else {
            0.0
        };
        transition_probabilities.insert(pair, prob);
    }

    transition_probabilities
}

fn main() {
    // Expect command line arguments:
    // 1) input FASTQ
    // 2) PWM output file
    // 3) MRF output file
    let args: Vec<String> = env::args().collect();
    if args.len() != 4 {
        eprintln!("Usage: {} <input_fastq> <pwm_output> <mrf_output>", args[0]);
        std::process::exit(1);
    }

    let input_fastq = &args[1];
    let pwm_output_path = &args[2];
    let mrf_output_path = &args[3];

    // Read sequences from FASTQ
    let mut seqs = Vec::new();
    let reader = fastq::Reader::from_file(input_fastq).expect("Could not open FASTQ file");
    for record in reader.records() {
        let rec = record.expect("Error reading record");
        seqs.push(rec.seq().to_vec());
    }

    // Convert sequence bytes to String
    let string_seqs: Vec<String> = seqs.iter()
                                       .map(|s| String::from_utf8_lossy(s).into_owned())
                                       .collect();

    // Build PWM
    let pwm = build_pwm(&string_seqs);

    // Build MRF (1st-order transitions)
    let mrf = build_mrf(&string_seqs);

    // Write PWM results
    let pwm_file = File::create(pwm_output_path).expect("Cannot create PWM output file");
    let mut pwm_writer = BufWriter::new(pwm_file);

    writeln!(pwm_writer, "Position Weight Matrix (probabilities)").unwrap();
    for (pos, position_map) in pwm.iter().enumerate() {
        writeln!(
            pwm_writer,
            "Position {}: A={:.3}, C={:.3}, G={:.3}, T={:.3}",
            pos,
            position_map[&'A'],
            position_map[&'C'],
            position_map[&'G'],
            position_map[&'T']
        ).unwrap();
    }

    // Write MRF results
    let mrf_file = File::create(mrf_output_path).expect("Cannot create MRF output file");
    let mut mrf_writer = BufWriter::new(mrf_file);

    writeln!(mrf_writer, "1st-order Markov Random Field (transition probabilities)").unwrap();
    for base1 in &['A', 'C', 'G', 'T'] {
        for base2 in &['A', 'C', 'G', 'T'] {
            let probability = mrf[&(*base1, *base2)];
            writeln!(mrf_writer, "{}->{}: {:.4}", base1, base2, probability).unwrap();
        }
    }

    println!("PWM and MRF computation completed successfully!");
}
