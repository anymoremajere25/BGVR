## Multi-Omics Integration and Nextflow Pipelines (Experiment_10_5)

This experiment contains a production-grade pipeline for multi-omics integration, combining high-resolution ChIP-seq/ATAC-seq peak detection and RNA-seq-based expression profiling with robust statistical analysis, all orchestrated via Nextflow. The pipeline includes:

---

###  Project Structure

```text
.
├── multiomics_calc_rs/                 # Rust tool to calculate coverage from BAM files
├── multiomics_main.rs/                      # Rust tool to perform multi-omics correlation and integration
├── multiomics_dataset_generator.py     # Python script to generate synthetic genomes and reads
├── main.nf                  # Nextflow pipeline to coordinate the workflow
├── data/                    # Input FASTQ/BAM files and annotations
└──Cargo.toml                
```

---

### Pipeline Overview

This pipeline performs:

1. **Coverage Calculation**
   Rust-based tool `calc_rs` calculates coverage over annotated gene regions using BAM files. Supports methods like `mean`, `total`, `density`, with optional normalization (RPKM-style), duplicate filtering, MAPQ filtering, and promoter extension.

2. **Fragment Shift Estimation**
   Estimation of optimal fragment shift for ChIP/ATAC-seq using cross-correlation (optionally FFT-accelerated). Produces shift estimates with confidence scores and signal-to-noise ratio.

3. **Peak Calling**
   A shift-aware Rust-based peak caller that identifies statistically significant enriched regions using a Poisson-based model with multiple testing correction.

4. **Expression Quantification**
   RNA-seq expression levels are loaded from quantification tools (e.g. Salmon, FeatureCounts) and analyzed across conditions.

5. **Multi-Omics Integration**
   Rust tool `multiomics_rs` integrates coverage and expression profiles using Pearson, Spearman, or Kendall correlation. It outputs:

   * Gene-level multi-omics integration scores
   * Correlation matrices and GraphML networks
   * Summary statistics and QC metrics

6. **Synthetic Dataset Generator**
   Python script generates synthetic genomes, GTF annotations, ChIP-seq reads, and RNA-seq expression matrices for benchmarking and validation.

---

### ⚙️ Technologies Used

| Language/Tool    | Role                                                                  |
| ---------------- | --------------------------------------------------------------------- |
| Rust             | High-performance computation for coverage, peak calling, and analysis |
| Nextflow         | Workflow orchestration with multi-step process control                |
| Python           | Simulation and data generation for testing                            |
| Samtools/Bowtie2 | Alignment and file preprocessing                                      |
| Polars           | High-speed tabular data processing in Rust                            |
| Rayon            | Thread-safe parallel processing                                       |
| Clap/Serde       | CLI parsing and config serialization                                  |

---

### 📝 Example Usage

#### 1. Run the Full Pipeline with Nextflow

```bash
nextflow run main.nf \
  --input_dir './data' \
  --output_dir './results' \
  --reference_genome './reference/genome.fa'
```

#### 2. Run Coverage Calculation Independently

```bash
cargo run --release --bin calc_rs -- \
  --bam_file sample.bam \
  --annotation genes.gtf \
  --output sample_coverage.csv \
  --method mean \
  --normalize_length \
  --extend_regions 2000 \
  --min_mapq 10 \
  --unique_only \
  --verbose
```

#### 3. Run Multi-Omics Integration

```bash
cargo run --release --bin multiomics_rs -- \
  --coverage coverage.csv \
  --expression expression.csv \
  --annotation genes.gtf \
  --correlation_method pearson \
  --output integration_results.json \
  --min_correlation 0.3 \
  --pvalue_threshold 0.05 \
  --output_matrix full_matrix.csv \
  --output_network network.graphml
```

---

## 📊 Output Summary

* `coverage.csv` – Per-gene coverage data
* `*_shift_estimate.json` – Fragment shift and correlation profile
* `*_peaks.bed` – Peak calls (with scores and significance)
* `integration_results.json` – Multi-omics correlation and scores
* `peak_calling_report.html` – HTML summary dashboard
* `merged_peaks.bed` – Aggregated peaks across samples
* `shift_comparison.txt` – Shift consistency report

---

### 📈 Visualization and Reporting

* Interactive HTML dashboards for QC and summary
* GraphML network for significant gene-gene co-regulation
* JSON reports for reproducibility and data mining

---

### 🧪 Simulated Test Data

To create a synthetic dataset for benchmarking:

```bash
python3 dataset_generator.py
```

Generates:

* Simulated `genome.fa`
* GTF annotations
* Synthetic ChIP-seq/ATAC-seq SAM files
* RNA expression matrices

---

### 🧠 Applications

* Chromatin accessibility and TF binding analysis
* Gene expression vs chromatin state association
* Multi-condition epigenomic profiling
* Benchmarking of peak callers and coverage metrics

---

### 📚 Citation and Attribution

If you use this pipeline in your work, please cite:

> Azkadina, M. (2025). Multi-Omics Integration and Nextflow Pipelines for Epigenomic Analysis. *Unpublished project report*.

---


