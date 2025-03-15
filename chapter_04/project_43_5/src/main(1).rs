use rayon::prelude::*;
use bio::io::fasta;
use std::sync::{Arc, Mutex};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::env;

// We'll use Serde for structured JSON output
use serde::Serialize;

#[derive(Debug, Serialize)]
struct MotifHit {
    seq_id: String,
    position: usize,
    score: f64,
}

/// A basic Position Weight Matrix structure for demonstration.
struct Pwm {
    // For each position, store log-odds or raw scores for A, C, G, T
    matrix: Vec<[f64; 4]>,
}

impl Pwm {
    /// Scores a DNA window against the PWM. Returns NEG_INFINITY for invalid bases.
    fn score_window(&self, window: &[u8]) -> f64 {
        let mut total_score = 0.0;
        for (pos, &base) in window.iter().enumerate() {
            let col_index = match base {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => { return f64::NEG_INFINITY; }
            };
            // Add the position-specific base score
            total_score += self.matrix[pos][col_index];
        }
        total_score
    }

    fn length(&self) -> usize {
        self.matrix.len()
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // 1) Parse command-line arguments
    //    cargo run --release <fasta_chunk> <output_json>
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: {} <fasta_chunk> <output_json>", args[0]);
        std::process::exit(1);
    }
    let fasta_path = &args[1];
    let output_json = &args[2];

    // 2) Define a sample PWM (6-bp motif). In real usage, read from config or a known motif model.
    let pwm = Pwm {
        matrix: vec![
            [1.2, -0.3, -0.3, -0.5],
            [0.5,  0.1, -0.2,  0.0],
            [1.1, -0.4, -0.6, -0.3],
            [0.7,  0.2, -0.1, -0.3],
            [0.9, -0.2, -0.3, -0.4],
            [0.3, -0.1,  0.2, -0.3],
        ]
    };

    let motif_len = pwm.length();

    // 3) Open the chunked FASTA
    let reader = fasta::Reader::from_file(fasta_path)?;

    // 4) Use a thread-safe vector to store hits
    let hits = Arc::new(Mutex::new(Vec::new()));

    // 5) Parallelize over FASTA records using rayon's par_bridge
    reader.records().par_bridge().for_each(|result| {
        let record = result.expect("Error reading FASTA record");
        let seq = record.seq();
        let seq_len = seq.len();

        // We'll keep local hits in a temporary vector to minimize lock contention
        let mut local_hits = Vec::new();
        for i in 0..(seq_len.saturating_sub(motif_len) + 1) {
            let window = &seq[i..(i + motif_len)];
            let score = pwm.score_window(window);
            // Arbitrarily define a threshold
            if score > 1.0 {
                local_hits.push(MotifHit {
                    seq_id: record.id().to_string(),
                    position: i,
                    score,
                });
            }
        }

        // Merge local hits into the global hits vector
        let mut global_hits = hits.lock().unwrap();
        global_hits.extend(local_hits);
    });

    // 6) Write partial hits to JSON
    let final_hits = hits.lock().unwrap();
    let out_file = File::create(output_json)?;
    let mut writer = BufWriter::new(out_file);
    for hit in final_hits.iter() {
        // Use Serde for structured JSON lines
        let json_line = serde_json::to_string(hit)
            .unwrap_or_else(|_| String::from("{\"error\": \"serialize\"}"));
        writer.write_all(json_line.as_bytes())?;
        writer.write_all(b"\n")?;
    }

    Ok(())
}
