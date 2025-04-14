### 6.6 Quality Control and Error Modeling
## experiment 6.6.1

The project successfully processes BAM files to compute quality control (QC) statistics using a Rust program (`coverage_tool`) and orchestrates the workflow with Nextflow. The Rust code, leveraging `rust-htslib` and `rayon`, analyzes a synthetic BAM file (`test_data/test.bam`, generated from `test.sam` with three reads on `chr1`), producing a JSON output (`qc_test.bam.json`) with coverage and mismatch counts (e.g., `{"coverage":3,"mismatches":3}`) for the region `chr1:1-300`. The Nextflow pipeline (`main.nf`) reads a sample list (`samples.txt`), executes the `collectQC` process to run `coverage_tool`, and stores the JSON output in a `work/` subdirectory, enabling scalable QC analysis for genomic data.

---

