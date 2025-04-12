## 6.5.  Advanced Data Structures for HTS Analysis

### Experiment_6_5_1

The process involves a genomic analysis pipeline using a Rust program (`genomic_interval_tree`) and Nextflow (`workflow.nf`) to compute and query genomic coverage from BAM-like files. Here's a detailed breakdown:
The process involves a genomic analysis pipeline using a Rust program (`genomic_interval_tree`) and Nextflow (`workflow.nf`) to compute and query genomic coverage from BAM-like files. Here's a detailed breakdown:

1. **Rust Program Overview (`main.rs`)**: The program has two main commands—`Coverage` and `Query`—built with `clap` for CLI parsing. It uses an interval tree data structure to manage genomic intervals efficiently.
   - **Coverage Command**: Reads a BAM-like file (e.g., `dummy1.bam.txt`, `dummy2.bam.txt`), parses start-end positions, and computes coverage by grouping overlapping reads into intervals with a `chunk_size` (default 50,000). Outputs a TSV file with `start`, `end`, and `coverage` (e.g., `coverage_dummy1.tsv`: `10 25 2`).
   - **Query Command**: Builds an interval tree from a coverage file, then queries it with a genomic range (e.g., `chr1:12-18`) to find overlapping intervals. Outputs results to a TSV (e.g., `query_result_chr1:12-18.tsv`).

2. **Nextflow Workflow (`workflow.nf`)**: Orchestrates the pipeline across multiple BAM files and queries.
   - **Input Files**: `bams.txt` lists BAM files (`dummy1.bam.txt`, `dummy2.bam.txt`), and `genome_intervals.txt` lists query intervals (`chr1:12-18`, `chr1:22-30`).
   - **Coverage Computation**: The `coverageComputation` process runs the Rust program’s `coverage` command on each BAM file, producing coverage TSVs (e.g., `coverage_dummy1.bam.txt.tsv` with `10 25 2`).
   - **Merge Coverage**: The `mergeCoverage` process concatenates all coverage TSVs into `merged_coverage.tsv`.
   - **Interval Query**: The `intervalQuery` process runs the Rust program’s `query` command on `merged_coverage.tsv` for each query interval, producing results like `query_result_chr1:12-18.tsv` (e.g., `10 25 2` for `chr1:12-18`).

3. **Execution and Results**:
   ![image](https://github.com/user-attachments/assets/ad06c92e-0778-4b37-80eb-a7b9d95da9da)

   - The screenshot shows Nextflow (version 24.10.5) running on a local executor with 5 tasks. Two BAM files are processed (`dummy1.bam.txt`, `dummy2.bam.txt`), each producing coverage TSVs (2 of 2 completed). These are merged (1 of 1 completed), and two queries (`chr1:12-18`, `chr1:22-30`) are executed (2 of 2 completed).
   - The result (`results.tsv`: `10 25 2`) indicates that the query `chr1:12-18` found an interval from position 10 to 25 with a coverage of 2, matching the merged coverage data.

In summary, the pipeline processes BAM files to compute coverage, merges results, and queries specific genomic intervals using an interval tree, all orchestrated by Nextflow for scalability and parallelism.
