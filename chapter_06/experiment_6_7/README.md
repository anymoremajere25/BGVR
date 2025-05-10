### 6.7. Integrative Analyses Using Rust-HTSlib

**Project Directory: experiment_6_7**

This section presents a Rust-based program that performs integrated analysis of read coverage and variant annotation. It showcases how developers can effectively utilize **rust-htslib** for reading BAM files, parse variants from BCF files, and enrich variant data with gene annotations from GFF files—all while leveraging Rust’s powerful concurrency model, logging infrastructure, and robust error handling. For more advanced tasks, crates like `ndarray` or `tch-rs` can be incorporated to support numerical computing or deep learning workflows. Rust’s strict type safety and memory model help prevent data corruption, a common challenge in large-scale genomics.

The program begins by parsing command-line arguments using the `clap` crate, which allows flexible specification of input files (BAM, BCF, and GFF). The function `load_gff_annotations` loads GFF annotation data into a `HashMap`. While this example uses a simple structure, real-world pipelines often use interval trees or suffix arrays for optimized querying.

After loading variant records from the BCF file, the code utilizes **Rayon** for parallel iteration. Each variant is processed concurrently by opening a new `IndexedReader` in the `process_variant` function to retrieve read coverage. Thanks to Rust's default immutability and thread safety, race conditions are avoided. Errors from individual variants are logged, but don’t crash the pipeline—enhancing robustness when dealing with noisy or incomplete datasets.

In production pipelines, integration with tools like **Nextflow** or **Snakemake** is common. Rust complements these with reliable performance and safe concurrency, suitable even for massive datasets. Libraries like `ndarray` or `tch-rs` enable advanced analysis and machine learning-based prioritization of variants. Rust’s guarantees ensure consistent behavior across multi-node clusters, handling terabytes of data efficiently.

---

### Nextflow Integration

The following **Nextflow DSL2** script extends the Rust pipeline into a modular, scalable workflow that can run on local machines, HPC clusters, or in the cloud. It divides the BCF file by chromosome, processes each segment with the BAM and GFF files using the Rust tool, and finally merges the results.

**1. splitBcf Process**

* **Input:** `cohort.bcf`
* **Action:** Indexes and splits the BCF by chromosome using `bcftools`.
* **Output:** Files like `split_chr1.bcf`.

**2. integrateData Process**

* **Inputs:**

  * Chromosome-specific BCF (`split_chr1.bcf`)
  * BAM file (`test1.bam`)
  * Annotation file (`annotations.gff`)
* **Action:**

  * Executes `rust_integrate_tool`.
  * Parses GFF into a lookup table.
  * Reads all variants and queries BAM coverage at each variant location.
  * Matches variants with gene annotations.
* **Output:** JSON files such as `integrated_test1.bam_split_chr1.bcf.json`, e.g.:

```
Variant: chr1:4 ref=A alt=T coverage=1 annotation=Some("ID=gene1;Name=GeneA")
```

**3. mergeIntegrations Process**

* **Input:** One or more JSON outputs from `integrateData`.
* **Action:** Merges them into a single output `final_integration.json`.

**Final Output:**

```json
[
  {
    "chrom": "chr1",
    "pos": 4,
    "ref_allele": "A",
    "alt_allele": "T",
    "coverage": 1,
    "annotation": "ID=gene1;Name=GeneA"
  }
]
```

---

### File Structure (excerpt)

```
experiment_6_7/
├── Cargo.toml
├── src/
│   ├── main.rs           # Rust source code
│   ├── main.nf           # Nextflow script
│   ├── test1.bam         # BAM file
│   ├── cohort.bcf        # BCF file
│   ├── annotations.gff   # GFF annotation file
│   ├── output.txt        # Output log
├── target/debug/
│   └── rust_integrate_tool.rar
├── work/.../
│   └── integrated_test1.bam_split_chr1.bcf.json
```

---

### Running the Workflow

**To run the Rust script manually in WSL:**

```bash
cargo run -- \
  --bam /mnt/c/.../test1.bam \
  --bcf /mnt/c/.../cohort.bcf \
  --gff /mnt/c/.../annotations.gff \
  | tee output.txt
```

**To execute the Nextflow pipeline:**

```bash
nextflow run main.nf \
  --bam_list 'bams.txt' \
  --bcf_file 'cohort.bcf' \
  --gff_file 'annotations.gff' \
  --bam_dir '.'
```

---

### ✅ Summary

* The Rust tool successfully annotated variants with gene info and computed read coverage.
* The Nextflow pipeline orchestrated data partitioning, parallel integration, and final result merging.
* Final output (`final_integration.json`) is ready for downstream tasks such as visualization or ML pipelines.
* This framework is reproducible, modular, and scalable for high-throughput genomics application.
