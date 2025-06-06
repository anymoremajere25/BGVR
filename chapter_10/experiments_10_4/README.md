## Experiment 10.4 - ChIP-Seq and ATAC-Seq Specific Algorithms

### Overview

This experiment showcases advanced fragment shift analysis and peak calling strategies specifically tailored for ChIP-Seq and ATAC-Seq datasets. Accurate estimation of fragment shifts using cross-correlation significantly improves peak localization. The components in this experiment integrate FFT acceleration, shift-aware peak calling, and high-performance Rust-based processing with a robust and scalable Nextflow pipeline.

### Objectives

* Perform strand cross-correlation analysis to estimate fragment shift sizes.
* Apply shift correction in peak calling to improve binding site localization.
* Validate results using simulated datasets with known shift properties.
* Ensure high-performance execution with multithreaded Rust code.
* Provide reproducible and interpretable outputs via an automated pipeline.

---

### Components

#### 1. `main.rs` (Rust Cross-Correlation Analyzer)

* Estimates fragment shift via strand cross-correlation (forward vs reverse reads).
* Supports FFT acceleration for up to 10x performance boost.
* Outputs quality metrics: NSC, RSC, estimated read length, SNR.
* Includes memory-efficient smoothing, binning, and sampling.
* JSON export of shift estimates with detailed coverage/correlation profiles.

#### 2. `peak_caller.rs` (Shift-Aware Peak Caller)

* Incorporates estimated shift into peak definition.
* Supports statistical validation with Poisson distribution.
* Outputs peaks in BED format with annotations and q-values.
* Filters based on peak width, significance, and coverage.

#### 3. `main.nf` (Nextflow Pipeline)

* Multi-sample support with consistent shift estimation.
* Quality-optimized alignment using `fastp` and duplicate removal.
* Automatically runs cross-correlation and peak calling.
* Generates HTML dashboard with shift profiles, metrics, and peak stats.

#### 4. `generate_shift_test_data.py` (Synthetic Dataset Generator)

* Simulates realistic ChIP-Seq data with specific fragment shifts (e.g., 100–250 bp).
* Embeds transcription factor binding site characteristics.
* Configurable GC content, enrichment ratio, fragment size.
* Generates gzipped FASTQ and BAM files for downstream testing.

#### 5. `shift_execution_guide.md`

* Step-by-step setup for running the full pipeline under WSL/Linux.
* Covers dataset generation, compilation, execution, and visualization.
* Includes troubleshooting for memory, BAM compatibility, and path issues.

---

### How to Run

#### Environment Setup

```bash
bash setup_environment.sh
```

#### Generate Synthetic Data

```bash
python3 generate_shift_test_data.py \
  --chromosomes chr1:1_000_000 chr2:1_000_000 \
  --fragment-shifts 100 150 200 250 \
  --output-dir test_data/
```

#### Run the Pipeline

```bash
nextflow run main.nf \
  --reads 'test_data/*.fastq.gz' \
  --genome test_data/ref.fa \
  --output results/
```

---

### Expected Output

* `shift_estimate.json`: Estimated shift, RSC, NSC, read length.
* `coverage_profiles.tsv`: Forward and reverse read coverage per bin.
* `correlation_profiles.tsv`: Raw and smoothed cross-correlation.
* `peaks.bed`: Shift-corrected significant peaks with annotations.
* `report.html`: Summary dashboard with shift plots and peak statistics.

---

### Key Metrics & Validation

* **Estimated Fragment Shift**: Close to ground truth (±10 bp).
* **RSC/NSC**: Validated against `phantompeakqualtools` benchmarks.
* **Peak Accuracy**: >90% overlap with known binding sites in test data.

---

### Technologies Used

* **Rust**: High-performance core processing
* **Rayon**: Parallel execution
* **RustFFT**: Accelerated Fourier transforms
* **Nextflow**: Workflow orchestration
* **FastQC/fastp**: Read quality control
* **SAMtools**: BAM processing

---

### Citation & References

* Zhang Y. et al., (2008). Model-based Analysis of ChIP-Seq (MACS).
* Landt S.G. et al., (2012). ChIP-seq guidelines and practices.
* phantompeakqualtools: [https://code.google.com/archive/p/phantompeakqualtools/](https://code.google.com/archive/p/phantompeakqualtools/)

---


