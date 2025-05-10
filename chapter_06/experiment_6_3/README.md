### 6.3. Variant Call Format (VCF/BCF) Handling

**Project: `experiment_6_3`**

This example demonstrates a production-ready approach to filtering BCF files in parallel using Rust. It combines `rust-htslib`, `rayon`, `clap`, `anyhow`, and `env_logger` for efficient data processing, error handling, and logging. The tool processes genomic data in configurable chunks to ensure memory efficiency, making it ideal for population-scale datasets in HPC or cloud environments.

#### Highlights:

* **Parallel Processing:** Each chunk of records is filtered concurrently using Rayon, effectively utilizing multicore systems.
* **Robust CLI:** Built with `clap` for flexible command-line argument parsing.
* **Safe & Efficient:** Rust's safety guarantees eliminate data races and memory errors.
* **Scalable:** Designed to handle massive variant datasets reliably across distributed systems.

#### Workflow Overview:

A corresponding **Nextflow** pipeline (`main.nf`) complements the Rust tool. It slices the input BCF file by genomic intervals using `bcftools`, then pipes each segment into the Rust-based `bcf_filter_tool`. Each chunk is handled as a separate task in a containerized environment, enabling high throughput and efficient resource use.

Filtered outputs per chunk (e.g., `filtered_chr1_1000-1200.bcf`) can later be merged with `bcftools concat` or custom merging steps.

#### Project Structure:

```
experiment_6_3/
├── Cargo.toml                   # Rust dependencies
├── src/
│   ├── main.rs                  # Rust filtering script
│   ├── main.nf                  # Nextflow pipeline
│   └── chunks_list.txt          # List of genomic regions
├── wgs_cohort.bcf               # Input BCF file
├── wgs_cohort.bcf.csi           # BCF index file
├── filtered_output.bcf          # Output from main.rs
└── target/debug/bcf_filter_tool.rar  # Compiled filter tool
```

#### How to Run:

**Standalone Rust Tool (full file filtering):**

```bash
cargo run -- --input wgs_cohort.bcf --output filtered_output.bcf --min-qual 30.0 --min-depth 10 --chunk-size 50000
```

* Filters all records with `QUAL ≥ 30` and average `DP ≥ 10`.
* Processes batches of 50,000 records in parallel.

**Nextflow Pipeline (region-based filtering):**

```bash
nextflow run main.nf
```

* Processes genomic intervals defined in `chunks_list.txt`.
* Each chunk runs independently in a parallel task using the Rust binary.

#### Output Summary:

| Execution              | Input         | Output                 | Parallelism       | Use Case                            |
| ---------------------- | ------------- | ---------------------- | ----------------- | ----------------------------------- |
| `cargo run`            | Full BCF file | `filtered_output.bcf`  | Batch-based       | Standalone, benchmark, full dataset |
| `nextflow run main.nf` | Chunked BCF   | `filtered_<chunk>.bcf` | Per genomic chunk | Scalable regional analysis          |

#### Key Benefits:

* **Modular and Reusable:** Works independently or as part of a workflow.
* **Highly Scalable:** Efficient use of resources across clusters or cloud.
* **Future-Ready:** Extendable for ML integration using `tch-rs`, `ndarray`, or `linfa`.



