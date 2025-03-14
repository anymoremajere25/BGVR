## 4.4. EPIGENOMIC DATA INTEGRATION AND ALGORITHMS
### Explanation of the Code and Output [project_44_1]
### Overview

This Rust program processes chromosome coverage data, applies optional smoothing, detects peaks based on a threshold, and writes the results to a file. The computations are parallelized using the Rayon library for efficiency.
### Code Breakdown
1. Struct Definition: ChromCoverage

struct ChromCoverage {
    chrom: String,
    data: Vec<f64>,
}

    Defines a structure to store coverage data for a chromosome.
    chrom: Name of the chromosome (e.g., "chr1", "chr2").
    data: A vector of floating-point numbers representing coverage values.

2. smooth_coverage Function (Optional Smoothing)

fn smooth_coverage(data: &[f64], window: usize) -> Vec<f64> {

    Applies a rolling average (window-based smoothing).
    Uses prefix sum for an efficient O(n) complexity.
    If window <= 1, it returns the original data.
    Computes the smoothed value at each position using the average of a local window.

3. local_peak_call Function (Detecting Peaks)

fn local_peak_call(data: &[f64], window: usize, threshold: f64) -> Vec<(usize, f64)> {

    Identifies peaks by computing local averages within a window.
    If the local average exceeds the threshold, the position is recorded as a peak.

Steps:

    Iterates through each index in data.
    Computes the local mean over a window.
    If the mean exceeds threshold, stores the (index, mean) pair.

4. call_peaks_and_smooth Function (Parallel Processing)

fn call_peaks_and_smooth(
    coverages: Vec<ChromCoverage>,
    window: usize,
    threshold: f64,
    do_smooth: bool,
) -> Vec<(String, usize, f64)> {

    Processes multiple chromosomes in parallel using Rayon.
    For each chromosome:
        Applies optional smoothing.
        Calls local_peak_call to detect peaks.
        Returns a list of tuples: (chromosome, position, peak_value).

Parallelization:

    into_par_iter() enables parallel processing.
    Each chromosome is processed independently.

5. main Function (Execution and File Output)

fn main() -> Result<(), Box<dyn std::error::Error>> {

    Defines sample chromosome coverage data:

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

    Sets:
        window = 3
        threshold = 3.0
        do_smooth = true (smoothing is applied)

    Calls call_peaks_and_smooth() to detect peaks.

    Writes output to partial_peaks.bed file.

### Output Explanation

The results are stored in partial_peaks.bed:

chr2	0	4.458
chr2	1	5.383
chr2	2	5.944
chr2	3	6.189
chr2	4	6.189
chr2	5	6.611
chr2	6	6.428
chr2	7	6.025

How the Output is Generated

    Only chr2 has peaks because its values are higher than the threshold after smoothing.
    chr1 is missing because none of its smoothed values exceed the threshold (3.0).
    The format:

    <chromosome>  <position>  <smoothed_peak_value>

### Key Takeaways

    Parallel Processing: Rayon optimizes computations across chromosomes.
    Efficient Smoothing: Uses prefix sum for O(n) complexity.
    Local Peak Detection: Peaks are identified based on local mean.
    File Output: Saves results in BED file format (common in genomics).

