# 🧬 experiment\_10\_1 — Fundamentals of Epigenomics

**Epigenomic Peak Calling with Rust, Python & Nextflow**

---

## 📖 Overview

This experiment introduces the fundamental techniques in epigenomic data analysis by implementing a complete, reproducible pipeline for peak detection in synthetic ChIP-Seq/ATAC-Seq data. We integrate a **Rust-based peak caller**, **Python data generation**, and **Nextflow orchestration** to simulate a realistic workflow. All dependencies are automatically installed and optimized for WSL/Linux.

---

## 🎯 Goals

* Generate synthetic FASTQ-like read profiles with known enriched regions
* Build a performant peak caller using Rust and Poisson statistics
* Construct a modular pipeline using Nextflow DSL2
* Automate environment setup for immediate reproducibility
* Output interpretable peak statistics and final HTML reports

---

## 📁 File Structure

| File                                               | Description                                                                                 |
| -------------------------------------------------- | ------------------------------------------------------------------------------------------- |
| [`Cargo.toml`](./Cargo.toml)                       | Rust project manifest with all dependencies                                                 |
| [`main.rs`](./src/main.rs)                         | Rust-based multithreaded peak caller with error handling and JSON/BED output                |
| [`main.nf`](./main.nf)                             | Nextflow DSL2 pipeline: QC, trimming, alignment, peak calling, reporting                    |
| [`generate_test_data.py`](./generate_test_data.py) | Python script to generate synthetic reads with simulated enrichment                         |
| [`setup_environment.sh`](./setup_environment.sh)   | WSL/Linux automation script for installing Rust, Python, bioinformatics tools, and Nextflow |
| [`run_commands.md`](./run_commands.md)             | Step-by-step execution instructions and usage examples                                      |

---

## ⚙️ Tools Used

* **Rust** (`rust-htslib`, `clap`, `rayon`, `serde`) for core peak detection
* **Python 3.9** for synthetic data generation & reporting
* **Nextflow DSL2** for pipeline orchestration
* **Bioinformatics stack**: `fastqc`, `fastp`, `bowtie2`, `samtools`, `bedtools`
* **WSL (Ubuntu)** for reproducible system setup

---

## 🔁 Pipeline Overview

```bash
FASTQ → QC → Trimming → Bowtie2 Alignment → BAM → Rust Peak Calling → JSON + BED + HTML
```

**Nextflow Stages:**

1. Input validation
2. Quality control (FastQC)
3. Adapter trimming (fastp)
4. Index building (bowtie2-build, samtools faidx)
5. Alignment (bowtie2 + samtools)
6. Peak calling (Rust binary with Poisson stats)
7. Merging and BED conversion
8. Final HTML report with peak statistics

---

## 🚀 How to Run

### ✅ 1. Setup the Environment

```bash
bash setup_environment.sh
source ~/.bashrc
```

This will:

* Install Rust, Python, conda, and all dependencies
* Configure aliases: `activate-biotools`, `epi-cd`
* Create `~/epigenomic_pipeline/` with proper directory structure

---

### 🔄 2. Run the Pipeline

```bash
cd ~/epigenomic_pipeline
./run_pipeline.sh
```

This will:

* Generate synthetic data in `epigenomic_test/`
* Compile the Rust peak caller
* Run the full Nextflow pipeline
* Output results to: `epigenomic_test/results/`

---

## 📤 Expected Outputs

| Output File                       | Description                                |
| --------------------------------- | ------------------------------------------ |
| `results/qc/*.html`               | FastQC reports                             |
| `results/alignments/*.bam`        | Sorted BAM alignments                      |
| `results/peaks/*.json`            | Peak regions in structured format          |
| `results/peaks/*.bed`             | Genome-browser-ready BED files             |
| `results/peaks/*.txt`             | Peak-level statistical summaries           |
| `results/merged/merged_peaks.bed` | Merged BED file across all samples         |
| `results/pipeline_report.html`    | Full HTML summary report (view in browser) |

---

## 📌 Parameter Highlights

The following parameters are tunable in `pipeline.config`:

| Parameter          | Default | Description                                   |
| ------------------ | ------- | --------------------------------------------- |
| `window_size`      | 200     | Sliding window size for scanning              |
| `min_coverage`     | 5.0     | Minimum read coverage to be considered a peak |
| `pvalue_threshold` | 0.05    | Significance threshold (Poisson model)        |
| `threads`          | 4       | Multithreading for each stage                 |

---

## 📚 What You’ll Learn

* Simulate epigenomic data and test peak callers against known truth
* Use Poisson statistics to evaluate signal over background
* Build and maintain reproducible pipelines using Nextflow
* Apply Rust for fast, memory-safe genome processing


