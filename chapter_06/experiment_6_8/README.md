### 6.8. Summary and Future Directions

**experiment_6_8**

This section presents a Rust-based command-line tool for calculating coverage from BAM files. The tool uses modern libraries including `clap` for command-line argument parsing, `anyhow` for error handling, and `rust-htslib` for working with BAM/CRAM files. Designed with scalability in mind, it fits seamlessly into HPC or cloud workflows managed by **Nextflow**. Containerization using Docker or Singularity ensures reproducibility and consistent environments across systems.

The tool reads a BAM file and genomic region provided by the user, computes the number of aligned reads in that region, and writes the result to a text file. It requires an indexed BAM file (.bai) for random region access. The combination of `anyhow` and `env_logger` supports robust error reporting and logging, which is particularly valuable in distributed environments.

When applied in parallelized workflows—such as splitting the genome into regions or processing multiple BAM files—the binary runs in ephemeral containers across compute nodes, each writing partial outputs. This approach reduces both memory use and total runtime.

For workflows incorporating machine learning, the `tch-rs` crate allows PyTorch integration, while `ndarray` supports vectorized numerical operations. Libraries like `polar`, `polars`, and `linfa` can be added for structured querying and traditional ML techniques. Despite increasing complexity, Rust’s type safety and concurrency features help maintain the reliability of genomic analysis pipelines.

The accompanying Nextflow workflow automates this process: each task handles a BAM-region pair in an isolated container, computes coverage using the Rust tool, and outputs results to individual text files. These partial results are then merged in a final step. This architecture supports high-throughput applications involving large sample sets or many genomic intervals, easily scaling across clusters or cloud services.

The `coverageCalc` process in Nextflow invokes the `rust_coverage_tool` (compiled from `main.rs`). Each container processes a unique combination of BAM file and genomic region, outputting results like `coverage_<bam>_<region>.txt`. Nextflow then collects and merges these outputs into `merged_coverage.txt`.

In production, additional features such as retries for transient failures, structured logging (e.g., for ELK stack integration), and version pinning for containers can be implemented to ensure reproducibility and traceability at scale.

#### File Structure

```
experiment_6_8/
  ├── Cargo.toml                 # Rust dependencies
  ├── src/
  │   ├── main.rs                # Rust tool source code
  │   ├── main.nf                # Nextflow workflow script
  │   ├── coverage_result.txt    # Example output
  │   ├── merged_coverage.txt    # Merged final output
  │   ├── bams.txt               # List of BAM files
  │   ├── regions.txt            # List of genomic regions
  │   ├── test.fa                # Reference FASTA
  │   ├── test1.bam              # Indexed BAM input
  │   ├── test1.bam.bai          # BAM index file
  │   ├── output.txt             # Additional output
  └── target/debug/
      └── rust_coverage_tool.rar # Compiled tool
```

#### How to Run

**Run `main.rs` in WSL:**

```bash
cargo run -- --bam test1.bam --region chr1:1-100 --out coverage_result.txt
```

**Run `main.nf` in WSL:**

```bash
nextflow run main.nf --bam_list bams.txt --region_list regions.txt --outdir .
```

#### Dependencies

```toml
anyhow = "1.0"
clap = { version = "4.4", features = ["derive"] }
log = "0.4"
env_logger = "0.11"
rust-htslib = "0.49.0"
```

---

### 🦀 `main.rs` – Rust Coverage Tool

**What it does:**

* Accepts three arguments: `--bam`, `--region`, and `--out`.
* Opens the BAM file and its index.
* Fetches all reads in the specified region.
* Counts aligned reads and writes the total to an output file.

**Example:**

```bash
cargo run -- --bam test1.bam --region chr1:1-100 --out coverage_result.txt
```

**Output:**

```
Computing coverage for region chr1:1-100 in BAM test1.bam
Coverage result: 1
```

Output file `coverage_result.txt` will contain:

```
1
```

---

### 🚀 `main.nf` – Nextflow Pipeline

**What it does:**

* Reads BAM paths from `bams.txt` and regions from `regions.txt`.
* Indexes BAMs using `samtools`.
* Creates BAM-region combinations.
* Runs `rust_coverage_tool` for each pair.
* Merges individual outputs into `merged_coverage.txt`.

**Example Inputs:**
`bams.txt`:

```
test1.bam
```

`regions.txt`:

```
chr1:1-100
chr1:200-300
```

→ This results in 2 runs:

* `coverage_test1.bam_chr1_1_100.txt` → contains 1
* `coverage_test1.bam_chr1_200_300.txt` → contains 0

Merged Output:

```text
# Merged coverage results
1
0
```

---

### ✅ Conclusion

* The Rust tool accurately computes coverage in specified genomic regions.
* The Nextflow pipeline successfully automates the process for multiple BAM files and regions.
* Output in `merged_coverage.txt` provides expected read counts per region.
* The modular and containerized design supports scaling in both HPC and cloud environments, ensuring reproducibility and efficient execution for large genomic datasets.


