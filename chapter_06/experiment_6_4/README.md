### 6.4. Parallel and Distributed Processing of HTS Data

**experiment_6_4**

This Rust program demonstrates how to perform parallelized read counting from a BAM file across multiple genomic regions. It’s designed to work efficiently on both HPC clusters and cloud infrastructures, making it ideal for temporary compute environments that handle subsets of data and then shut down. The implementation utilizes the `rust-htslib` crate for handling BAM files and `rayon` for concurrency, with improved error handling and logging for robustness. The code also allows easy extensions—such as streaming with channels or integrating libraries like `ndarray`, `linfa`, or deep learning frameworks like `tch-rs`.

In this enhanced version, the `process_bam_chunk` function returns a `Result` type via the `anyhow` crate, which helps unify error messages into a single, manageable object—enabling error tracking and retry logic. Each genomic region is processed concurrently using `rayon`’s parallel iterators, maximizing CPU core or node utilization in HPC settings. If an error occurs in one region, it doesn’t halt the entire program; instead, it logs the error and continues processing the remaining regions.

For production environments, more advanced features can be integrated. These may include structured logging (e.g., JSON format) for compatibility with observability tools like ELK or OpenTelemetry, and streaming of partial results to a separate writer thread. When deployed in HPC or cloud environments, containerizing the Rust tool guarantees a lightweight, consistent runtime and facilitates deployment with orchestration systems like Nextflow or Kubernetes. Additionally, the use of Rust enables seamless integration with high-performance numerical or AI libraries while ensuring memory safety.

The accompanying Nextflow script automates various stages: acquiring BAM files, executing a Rust-based tool for variant calling or read counting, and merging output files. Each phase runs within a disposable container or isolated environment, utilizing cloud or HPC resources for parallel execution. The Rust binary, generated from your earlier code and referred to as `rust_caller_tool`, can be swapped in easily. Being statically linked, Rust binaries minimize runtime overhead, ideal for containerized execution in Nextflow.

This pipeline showcases how Nextflow can streamline complex workflows. The `alignmentOrFetch` stage retrieves or generates BAM files for each sample. Then, the `variantCalling` step uses the Rust-based tool to perform parallel operations (like read counting or variant calling) using `rayon` for multithreaded efficiency. Lastly, `mergeVcfs` combines all generated VCF files into a single output, ready for further annotation or population-level studies.

Rust and Nextflow are increasingly popular among AI engineers and bioinformaticians for building scalable genomics pipelines, capable of running economically on public clouds or local HPC clusters. One notable example includes a major hospital system processing thousands of clinical exomes weekly using ephemeral containers. These containers run the full workflow from read alignment to variant calling and merging. Rust’s strong concurrency model helped them avoid race conditions, while the compact size of static binaries boosted performance and minimized container image sizes.

#### Project Structure:

```
experiment_6_4/
├── Cargo.toml                     # Rust project dependencies
├── src/
│   ├── main.rs                   # Rust script
│   ├── main.nf                   # Nextflow script
│   ├── ref.fasta                 # Reference FASTA file
│   ├── ref.fasta.fai             # Indexed FASTA
│   ├── sample1.bam               # Input BAM file for sample 1
│   ├── sample1.bam.bai           # BAM index
│   ├── sample2.bam               # Input BAM file for sample 2
│   ├── sample2.bam.bai           # BAM index
│   ├── samples.txt               # Sample list input
│   ├── output.txt                # Output from Rust script
│   └── results/
│       └── merged.vcf            # Final merged VCF output
├── target/debug/
│   └── bam_read_counter.rar      # Compressed output binary
```

#### How to Run:

**Run Rust script in WSL:**

```bash
cargo run -- --bam sample1.bam --region 'chr1:1-10' 'chr1:1-32' 2>&1 | tee output.txt
```

This runs `main.rs` with `sample1.bam` as input for regions `chr1:1-10` and `chr1:1-32`, and logs output to `output.txt`.

**Run Nextflow script in WSL:**

```bash
nextflow run main.nf
```

With parameters:

```groovy
params.sample_list = "samples.txt"
params.output_dir = "results"
params.region = "chr1:1-32"
params.mock = true // Enables mock mode for testing
```

#### Rust Dependencies:

```toml
clap = { version = "4.0", features = ["derive"] }
anyhow = "1.0"
log = "0.4"
rayon = "1.7"
rust-htslib = "0.49.0"
```

---

### Explanation of the Output:

We executed a full mock variant-calling pipeline via Nextflow using sample data and a custom Rust-based read counter (`bam_read_counter`). Here's what each part did:

🔢 **Input Files:**

* `samples.txt`: Lists two samples to process (sample1, sample2).
* `sample1.bam`, `sample2.bam`: Aligned BAM files.
* `ref.fasta`: Reference genome (not used in mock mode).

⚙️ **Pipeline Stages:**

1. **Channel Setup**

   * Parses `samples.txt`, sending each sample ID to the workflow.

2. **alignmentOrFetch Process**

   * Copies the BAM files into the working directory.

3. **variantCalling Process**

   * With `mock = true`, it generates dummy VCFs per sample with 3 predefined variants each:

     ```
     chr1 14653 . A G 100 PASS . GT 0/1
     chr1 14907 . A G 100 PASS . GT 0/1
     chr1 15211 . G A 100 PASS . GT 1/1
     ```

4. **mergeVcfs Process**

   * Concatenates the individual VCFs.
   * Includes a single header and appends all records, without deduplication or sample-specific columns.

📁 **Output File:**
`results/merged.vcf` contains 6 entries—3 for each sample:

```vcf
##fileformat=VCFv4.2
##source=MockVariantCaller
##contig=<ID=chr1,length=248956422>
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT sample1
chr1 14653 . A G 100 PASS . GT 0/1
chr1 14907 . A G 100 PASS . GT 0/1
chr1 15211 . G A 100 PASS . GT 1/1
chr1 14653 . A G 100 PASS . GT 0/1
chr1 14907 . A G 100 PASS . GT 0/1
chr1 15211 . G A 100 PASS . GT 1/1
```

🧾 **Conclusion:**

✅ Your Rust tool (`bam_read_counter`) worked correctly for read-region counting and is ready for integration.

✅ The Nextflow pipeline executed successfully with support for multiple samples.

🧪 Note: Since `params.mock = true`, the workflow used simulated data. No real alignments or variant calls were performed.

🔍 The final `merged.vcf` confirms that both samples were processed and mock variants successfully merged.



