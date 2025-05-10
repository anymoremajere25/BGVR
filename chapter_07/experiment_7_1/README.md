### 7.1. Foundational Data Structures in Rust Noodles

**experiment_7_1**

This experiment showcases the integration of Rust with the `noodles-bam` and `noodles-vcf` crates for basic genomic data processing—specifically, coverage calculation and variant file recognition. The Rust code is optimized for high-performance computing (HPC) by leveraging parallel execution with Rayon’s `.par_iter()` to process multiple genomic intervals concurrently.

Although the example uses a simple line-counting approach for coverage, the design is extensible. For instance, advanced implementations can maintain coverage in-memory using data structures like segment trees, enabling more sophisticated downstream analyses.

For numerical or AI-based processing, additional crates such as `ndarray` (for numerical arrays) or `tch-rs` (for PyTorch integration) can be used. This flexibility makes Rust suitable for both basic genomics and advanced bioinformatics pipelines.

The accompanying **Nextflow pipeline** demonstrates the orchestration of these Rust binaries in ephemeral container tasks. Each container processes a specific genomic region, computes partial results, writes outputs, and exits—ensuring efficient resource utilization and scalability in HPC or cloud environments.

The pipeline runs two Rust tools:

* `rust_noodles_coverage` for calculating coverage
* `rust_noodles_variant` for parsing variants (not fully implemented yet)

This containerized and ephemeral approach suits AI-driven genomics projects, where large datasets are split across multiple compute nodes. The use of fixed container images guarantees reproducibility, while Nextflow handles job dispatch, error management, and output merging.

---

**File Structure:**

```
experiment_7_1/
├── Cargo.toml                 # Rust dependencies
├── src/
│   ├── main.rs                # Rust script for coverage calculation
│   ├── main.nf                # Nextflow pipeline
│   ├── bams.txt               # List of sorted BAM files
│   ├── cohort.vcf             # Sample VCF file
│   ├── regions.txt            # Genomic regions list
│   ├── sample1.sam / .bam / .bai
│   ├── sample2.sam / .bam / .bai
│   ├── variants.vcf           # VCF file with variant calls
│   └── output.txt             # Output from Rust tool
├── target/debug/
│   └── rust_noodles_tool.rar  # Compressed binary output
└── work/
    └── coverage_output.txt    # Result from pipeline execution
```

---

**How to Run:**

In **WSL (Windows Subsystem for Linux)**:

1. Run the Rust script directly:

```bash
cargo run -- --vcf-file cohort.vcf --bam-files sample1.sorted.bam,sample2.sorted.bam | tee output.txt
```

This executes `main.rs` using `cohort.vcf` and two BAM files as input, and writes the output to `coverage_output.txt`.

2. Run the Nextflow pipeline:

```bash
nextflow run main.nf
```

Or, with explicit parameters:

```bash
nextflow run main.nf --bam_list "bams.txt" --vcf_file "cohort.vcf" --rust_bin "/mnt/c/Users/ragon/BGVR/chapter_07/experiment_7_1/target/debug"
```

---

**Rust Tool Overview (`main.rs`):**

This CLI tool performs basic BAM parsing using `noodles`. Key features:

* Uses `clap` for CLI argument parsing
* Initializes logging via `env_logger`
* Reads BAM headers and prints up to 2 reference sequences
* Processes BAM records in parallel using Rayon
* Prints limited details (e.g., position, MAPQ, CIGAR of first 5 records)
* Logs total records processed per BAM file

**Note:** The VCF file is currently acknowledged but not parsed—this feature can be enabled by adding `"vcf"` to the `noodles` crate features in `Cargo.toml`.

---

**Nextflow Pipeline (`main.nf`) Overview:**

* **Inputs:**

  * `bams.txt` (converted to comma-separated list)
  * `cohort.vcf`
  * Rust binary path

* **Main Process:**

  * Converts the BAM list to CLI format:

    ```bash
    BAM_FILES=$(cat bams.txt | tr '\n' ',' | sed 's/,$//')
    ```
  * Runs the compiled binary:

    ```bash
    rust_noodles_coverage --vcf-file cohort.vcf --bam-files "sample1.sorted.bam,sample2.sorted.bam"
    ```
  * Output is written to `coverage_output.txt`

---

**Sample Output:**

```
Starting BAM processing application...
VCF file noted: cohort.vcf (VCF processing not enabled)
Processing 2 BAM files...
Processing BAM file: sample1.sorted.bam
Processing BAM file: sample2.sorted.bam
Processing completed successfully
```

This minimal output is consistent with the `println!` statements in `main.rs`.

---

**Conclusion:**

✅ Rust tool successfully integrated with a Nextflow pipeline
✅ Input handling, parallelism, and logging are functional
⚠️ VCF support is acknowledged but not yet implemented
📁 Output is generated inside Nextflow’s working directory structure

