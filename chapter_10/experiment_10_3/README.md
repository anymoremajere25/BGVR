

## Experiment 10\_3 — Signal Normalization and Quality Control

This experiment implements a **comprehensive quality control (QC) and signal normalization framework** for ChIP-seq and ATAC-seq data using a hybrid toolchain built in **Rust**, **Python**, and **Nextflow**. It provides reproducible, scalable, and interpretable QC metrics crucial for downstream peak calling, normalization, and integration tasks.

---

### 📁 File Structure

```
experiment_10_3/
├── main.rs                    # Rust-based QC analyzer with FRiP, SNR, and normalization
├── peak_caller.rs             # Lightweight peak caller for pre-QC region generation
├── main.nf                    # Nextflow pipeline integrating QC, alignment, and reporting
├── generate_qc_test_data.py   # Python-based synthetic QC dataset generator
├── qc_execution_guide.md      # Setup, run instructions, and troubleshooting guide
└── output/                    # Results (BigWig, QC scores, plots, reports)
```

---

### 🎯 Goals

* Measure **quality metrics** such as:

  * **FRiP (Fraction of Reads in Peaks)**
  * **SNR (Signal-to-Noise Ratio)**
  * **Library complexity**, **duplication rate**
  * **Fragment size distribution**
* Normalize signals using:

  * **RPM (Reads Per Million)**
  * **RPKM (Reads Per Kilobase per Million)**
  * **TPM (Transcripts Per Million)**
* Generate **BigWig tracks** for visualization
* Identify high/low quality samples via thresholded pass/fail rules
* Visualize QC metrics across multiple samples

---

### 🔧 Components

### 1. 🔬 **Robust QC Analyzer in Rust (`main.rs`)**

* **Multi-threaded analysis** with Rayon for fast per-sample QC.
* Accepts BAM, BED, and optional blacklist/background regions.
* Computes:

  * FRiP, SNR, duplicate rate
  * Read complexity metrics
  * Signal normalization (RPM, RPKM, TPM)
* Output formats: `JSON`, `CSV`, `summary.txt`

### 2. ⚡ **Minimal Peak Caller (`peak_caller.rs`)**

* Sliding window-based region detector
* Ideal for generating approximate peaks for QC when true regions are unknown
* Outputs BED format compatible with the QC module

### 3. 🔁 **Nextflow Pipeline (`main.nf`)**

* Full QC pipeline: FASTQ → Trim → Align → QC → Normalize → Report
* Optional region generation (or uses known peak/background/blacklist files)
* **BigWig** generation for IGV or UCSC Genome Browser
* **HTML dashboard** with QC summary and CSV exports
* Automatic **pass/fail validation** based on predefined thresholds

### 4. 🧬 **QC Dataset Generator (`generate_qc_test_data.py`)**

* Generates test datasets in three quality tiers:

  * **High quality**: high FRiP, low duplication
  * **Medium**: moderate FRiP and SNR
  * **Low quality**: poor enrichment, high noise/duplication
* Simulates fragment distributions, GC bias, and regulatory enrichment
* Produces:

  * BAM/FASTQ reads
  * Known peaks, background, and blacklist regions
  * Summary metadata for validation

### 5. 📘 **Setup & Execution Guide (`qc_execution_guide.md`)**

* Environment setup (WSL + conda + Rust)
* Component-level testing instructions
* Performance tuning tips
* Guidelines for real dataset deployment

---

## 📤 Expected Output

| File                    | Description                                      |
| ----------------------- | ------------------------------------------------ |
| `qc_scores.json`        | Per-sample metrics: FRiP, SNR, duplication, etc. |
| `signal_normalized.csv` | RPM/RPKM/TPM values per region                   |
| `*.bw`                  | BigWig tracks for genome visualization           |
| `summary_report.html`   | Interactive HTML report (multi-sample view)      |
| `failures.csv`          | List of failed samples with reasons              |
| `fragment_dist.png`     | Fragment size histogram per sample               |
| `region_stats.txt`      | Summary of background vs enriched regions        |

---

## 📈 Key QC Metrics

| Metric                 | Description                                       |
| ---------------------- | ------------------------------------------------- |
| **FRiP**               | Fraction of total reads in peak regions           |
| **SNR**                | Signal-to-noise ratio (peak vs background)        |
| **Duplication Rate**   | Redundancy of mapped reads                        |
| **Library Complexity** | Unique vs total reads                             |
| **Peak Coverage**      | Coverage in known/putative enriched regions       |
| **Blacklist Overlap**  | Reads in low-confidence or artifact-prone regions |

---

## ⚙️ Running the Pipeline

### Step 1: Generate QC Datasets

```bash
python generate_qc_test_data.py --output qc_test
```

### Step 2: Build Rust Tools

```bash
cd rust_qc
cargo build --release
```

### Step 3: Run Nextflow

```bash
nextflow run main.nf \
  --input_dir qc_test \
  --output_dir qc_test/results \
  --threads 4
```

---

## 🧠 Interpretation Tips

* Use the HTML dashboard to explore per-sample FRiP and SNR.
* Compare BigWig tracks across quality tiers.
* Identify outliers via clustering of normalized signal matrices.
* Use fail/pass QC gating to exclude poor quality samples in downstream peak calling (10.4).

---

## ✅ Requirements

* **Rust** (1.70+)
* **Python 3.8+**
* **Nextflow 22.10+**
* `samtools`, `bedtools`, `fastqc`, `fastp`, `deeptools`
* Optional: `IGV`, `UCSC Genome Browser` for BigWig viewing

---

## 🔍 Real-World Use Case

This QC infrastructure is suitable for:

* **Large-scale ChIP-seq/ATAC-seq projects** across tissues or timepoints
* **Pre-filtering** poor quality samples before expensive downstream analysis
* **Validating enrichment** around biological features (e.g., promoters)
* **Batch comparison** and reporting in consortia projects

---


