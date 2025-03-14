use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};

/// Stores coverage data for a given chromosome.
struct ChromCoverage {
    chrom: String,
    data: Vec<f64>,
}

/// Optionally smooth coverage using a rolling-window average.
/// Uses a prefix-sum approach for O(n) complexity.
fn smooth_coverage(data: &[f64], window: usize) -> Vec<f64> {
    if window <= 1 {
        return data.to_vec();
    }

    let mut prefix_sum = Vec::with_capacity(data.len() + 1);
    prefix_sum.push(0.0);

    for &val in data {
        prefix_sum.push(prefix_sum.last().unwrap() + val);
    }

    let half_window = window / 2;
    let n = data.len();
    let mut smoothed = vec![0.0; n];

    for i in 0..n {
        let start = i.saturating_sub(half_window);
        let end = (i + half_window).min(n - 1);
        let sum_slice = prefix_sum[end + 1] - prefix_sum[start];
        let length = (end - start + 1) as f64;
        smoothed[i] = sum_slice / length;
    }
    smoothed
}

/// A minimal local peak-calling function.
fn local_peak_call(data: &[f64], window: usize, threshold: f64) -> Vec<(usize, f64)> {
    let half_window = window / 2;
    let n = data.len();
    let mut peaks = Vec::new();

    for i in 0..n {
        let start = i.saturating_sub(half_window);
        let end = (i + half_window).min(n - 1);
        let segment_mean = data[start..=end].iter().sum::<f64>() / (end - start + 1) as f64;
        if segment_mean > threshold {
            peaks.push((i, segment_mean));
        }
    }
    peaks
}

/// High-level function to call peaks across multiple chromosomes in parallel.
/// 1) Optionally smooth coverage.
/// 2) For each chromosome, call peaks using a local function.
/// 3) Collect results in a single flattened vector (chrom, position, peak_value).
fn call_peaks_and_smooth(
    coverages: Vec<ChromCoverage>,
    window: usize,
    threshold: f64,
    do_smooth: bool,
) -> Vec<(String, usize, f64)> {
    coverages
        .into_par_iter()
        .flat_map(|cov| {
            let data = if do_smooth {
                smooth_coverage(&cov.data, window)
            } else {
                cov.data.clone()
            };

            local_peak_call(&data, window, threshold)
                .into_iter()
                .map(move |(pos, val)| (cov.chrom.clone(), pos, val))
                .collect::<Vec<_>>()
        })
        .collect()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let coverage_data = vec![
        ChromCoverage {
            chrom: "chr1".to_string(),
            data: vec![0.0, 2.5, 5.5, 2.2, 0.9, 4.1, 3.5, 0.7],
        },
        ChromCoverage {
            chrom: "chr2".to_string(),
            data: vec![0.0, 7.5, 8.0, 6.2, 2.1, 9.4, 10.2, 0.5],
        },
    ];

    let window = 3;
    let threshold = 3.0;
    let do_smooth = true;

    let all_peaks = call_peaks_and_smooth(coverage_data, window, threshold, do_smooth);

    let file = File::create("partial_peaks.bed")?;
    let mut writer = BufWriter::new(file);

    for (chrom, pos, val) in all_peaks {
        let line = format!("{}\t{}\t{:.3}\n", chrom, pos, val);
        writer.write_all(line.as_bytes())?;
    }

    Ok(())
}
