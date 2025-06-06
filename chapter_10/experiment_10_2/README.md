# 🧬 Experiment 10\_2: Data Structures for Peak Calling

This experiment demonstrates a **complete interval-based peak calling system** implemented using a combination of **Rust**, **Nextflow**, and a **realistic synthetic dataset generator**. It is designed to simulate biologically plausible ChIP-seq data and perform statistically rigorous peak detection using advanced data structures like interval trees.

---

## 📁 Project Structure

```
experiment_10_2/
├── Cargo.toml                     # Rust crate definitions and dependencies
├── main.rs                        # Rust implementation for peak calling using interval trees
├── main.nf                        # Nextflow workflow for preprocessing, alignment, and calling
├── interval_execution_guide.md   # Step-by-step usage guide
├── generate_interval_test_data.py # Realistic synthetic dataset generator
└── README/                       # Output results from Nextflow pipeline (created after execution)
```

---

## 🚀 Key Features and Enhancements

### 1. Rust-Based Peak Caller (`main.rs`)

* **Efficient interval tree** for O(log n) lookup during peak matching.
* **Statistical modeling** using Poisson distribution with multiple testing correction.
* **Fragment shifting and extension** to mimic true ChIP-seq profiles.
* **Progress bar and logging** via `indicatif` and `log`.
* **Configurable BED input**, support for genome-wide or region-based calling.

### 2. Nextflow Pipeline (`main.nf`)

* **FastQC + fastp** for raw read QC and trimming.
* **Alignment with BWA**, duplicate marking, and filtering.
* **Interval management** with BED file fallback to dynamic generation.
* **Comprehensive reporting** (HTML summary, plots, stats).
* **Multi-sample support** and **cross-sample peak overlap** analysis.

### 3. Synthetic Dataset Generator (`generate_interval_test_data.py`)

* Generates:

  * A realistic reference genome with GC bias, repeat elements, CpG islands.
  * Biologically plausible **gene**, **promoter**, and **enhancer** intervals.
  * Simulated **ChIP-seq reads** enriched in regulatory regions.
  * True peak BED file for validation.

### 4. Execution Guide (`interval_execution_guide.md`)

* Instructions for:

  * Environment setup (WSL/Ubuntu).
  * Generating synthetic data.
  * Running pipeline end-to-end.
  * Interpreting output and validating performance.

---

## 📊 Expected Output

After executing the full pipeline, you will obtain:

| Output Type              | Description                            |
| ------------------------ | -------------------------------------- |
| `*.fastq.gz`             | Synthetic ChIP-seq reads per sample    |
| `genome.fa`              | Realistic reference genome             |
| `*.bed`                  | Gene, promoter, and enhancer intervals |
| `true_peaks.bed`         | Known ground truth peaks               |
| `peaks.json / peaks.bed` | Detected peaks from Rust binary        |
| `nextflow_report.html`   | QC and pipeline summary report         |
| `nextflow_trace.txt`     | Detailed execution trace               |
| `nextflow_timeline.html` | Resource usage timeline                |

You can compare detected peaks with `true_peaks.bed` to calculate:

* **Sensitivity (Recall)**
* **Specificity**
* **FDR / q-values**

---

## ⚙️ Quick Start Guide

### 1. Generate Synthetic Data

```bash
python generate_interval_test_data.py --output-dir ./interval_test
```

### 2. Build Rust Peak Caller

```bash
cd interval_test
cargo build --release
```

### 3. Run Pipeline

```bash
./run_pipeline.sh
```

OR run manually:

```bash
nextflow run ../main.nf -c nextflow.config
```

---

## 📈 Validation Strategy

The experiment enables quantitative evaluation using:

* Overlap between called peaks and known `true_peaks.bed`.
* Positional enrichment plots around TSS and enhancers.
* Peak count, width distribution, and fold-enrichment metrics.

---

## 🔍 Advanced Capabilities

| Capability                 | Description                              |
| -------------------------- | ---------------------------------------- |
| Interval trees             | Speed up genomic region queries          |
| Streaming I/O              | Reduces memory usage on large datasets   |
| Multi-threading (Rayon)    | Accelerates peak calling                 |
| Configurable via CLI/JSON  | Supports custom interval sets and tuning |
| BED or genome-wide support | Flexible input                           |
| Full HTML reports          | For easy visualization and debugging     |

---

## ✅ Requirements

* Rust (1.70+)
* Python 3.7+
* Nextflow
* BWA, Samtools, FastQC, fastp
* Conda environment (optional: `biotools`)

---

## 📚 References

* Interval Tree concept: [Cormen et al. "Introduction to Algorithms"](https://mitpress.mit.edu/9780262046305/)
* Poisson model for ChIP-seq: MACS2 framework
* Synthetic data generation: adapted from ENCODE and SimGenome approaches

---



