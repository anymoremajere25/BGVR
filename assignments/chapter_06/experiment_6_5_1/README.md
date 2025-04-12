### 6.5. Advanced Data Structures for HTS Analysis
#### Experiment_6_5_1

Based on the corrected codes we used (`main.nf`, `main.rs` for `rust_interval_query_tool`, and `coverage_tool/main.rs` for `rust_coverage_tool`). The pipeline simulates a bioinformatics workflow to analyze genomic coverage from BAM files, and I'll break down what it does step-by-step, what each part contributes, and what the output files mean, all tailored to your setup using WSL and Cursor. This will help the reader understand how the code works together to produce the results:

```
executor >  local (5)
[e9/6555a1] coverageComputation (mock_sample1.bam) | 2 of 2 ✔
[a8/5660e9] mergeCoverage                          | 1 of 1 ✔
[5c/f4e1a8] intervalQuery (chr1:1000-2000)         | 2 of 2 ✔
```

```

### Overview of the Pipeline
The pipeline mimics a real bioinformatics task: analyzing DNA sequencing data to find how much coverage (number of reads) exists in specific genomic regions. It uses:
- **Mock BAM files** (listed in `bams.txt`): Pretend sequencing data (`mock_sample1.bam`, `mock_sample2.bam`).
- **Genomic intervals** (in `genome_intervals.txt`): Regions to query (e.g., `chr1:1000-2000`).
- **Rust tools**:
  - `rust_coverage_tool`: Simulates calculating coverage from BAMs.
  - `rust_interval_query_tool`: Finds which coverage regions overlap the queried intervals.
- **Nextflow** (`main.nf`): Manages the workflow, running tasks in order or parallel as needed.

The pipeline has three main steps:
1. **Compute Coverage**: Process each BAM file to get coverage data.
2. **Merge Coverage**: Combine coverage from all BAMs into one file.
3. **Query Intervals**: Check which coverage regions overlap specific genomic intervals.

The output files show the coverage for each BAM and the overlaps for each queried region.

---

### Step-by-Step Process Explanation

#### 1. **Input Files**
The pipeline starts with two input files you created:
- **`bams.txt`**:
  ```
  mock_sample1.bam
  mock_sample2.bam
  ```
  - Lists two fake BAM files (sequencing data).
  - Each line is a file name that the pipeline will process.
  - Created with:
    ```bash
    echo "mock_sample1.bam" > bams.txt
    echo "mock_sample2.bam" >> bams.txt
    ```

- **`genome_intervals.txt`**:
  ```
  chr1:1000-2000
  chr1:1500-2500
  ```
  - Lists two genomic regions to query (chromosome 1, positions 1000-2000 and 1500-2500).
  - Format: `chr:start-end`.
  - Created with:
    ```bash
    echo "chr1:1000-2000" > genome_intervals.txt
    echo "chr1:1500-2500" >> genome_intervals.txt
    ```

- **How They’re Used**:
  - `main.nf` reads these files to create channels (like lists Nextflow processes):
    ```groovy
    bamChannel = Channel.fromPath(params.bam_list).splitText().map { it.trim() }
    queriesChannel = Channel.fromPath(params.ref_intervals).splitText().map { it.trim() }
    ```
  - `bamChannel`: Contains `mock_sample1.bam`, `mock_sample2.bam`.
  - `queriesChannel`: Contains `chr1:1000-2000`, `chr1:1500-2500`.
  - The `.map { it.trim() }` removes extra newlines to keep data clean.

#### 2. **Coverage Computation (Process: `coverageComputation`)**
- **What It Does**:
  - Takes each BAM file from `bamChannel` (e.g., `mock_sample1.bam`).
  - Runs `rust_coverage_tool` to simulate calculating coverage (how many reads cover genomic positions).
  - Produces a TSV file with coverage intervals for each BAM.
- **Code**:
  ```groovy
  process coverageComputation {
      tag "${bam_file}"
      input:
      val(bam_file)
      output:
      path "coverage_${bam_file}.tsv"
      script:
      """
      rust_coverage_tool \
          --bam "${bam_file}" \
          --out "coverage_${bam_file}.tsv" \
          --chunk-size ${params.parallel_chunck_size}
      """
  }
  ```
  - `tag "${bam_file}"`: Labels each task (e.g., `mock_sample1.bam`) for tracking.
  - `val(bam_file)`: Gets one BAM name from `bamChannel`.
  - `path "coverage_${bam_file}.tsv"`: Names the output (e.g., `coverage_mock_sample1.bam.tsv`).
  - Runs:
    ```
    rust_coverage_tool --bam "mock_sample1.bam" --out "coverage_mock_sample1.bam.tsv" --chunk-size 50000
    ```

- **Rust Tool (`rust_coverage_tool`)**:
  - Defined in `coverage_tool/src/main.rs`:
    ```rust
    #[derive(Parser, Debug)]
    struct Cli {
        #[arg(long)]
        bam: String,
        #[arg(long)]
        out: String,
        #[arg(long)]
        chunk_size: u64,
    }
    fn main() -> Result<()> {
        let cli = Cli::parse();
        let mut file = File::create(&cli.out)?;
        writeln!(file, "chr1\t1000\t2000\t50")?;
        writeln!(file, "chr1\t1500\t2500\t60")?;
        println!("Mock coverage written to {} for BAM {}", cli.out, cli.bam);
        Ok(())
    }
    ```
  - **What It Does**:
    - Takes `--bam` (e.g., `mock_sample1.bam`), `--out` (e.g., `coverage_mock_sample1.bam.tsv`), and `--chunk-size` (50000, unused in mock).
    - Writes fake coverage data to the output file:
      ```
      chr1    1000    2000    50
      chr1    1500    2500    60
      ```
    - Format: `chromosome\tstart\tend\tcoverage`.
    - Simulates real coverage (e.g., 50 reads from positions 1000-2000 on chr1).

- **Output**:
  - Two files (one per BAM):
    - `output/work/*/coverage_mock_sample1.bam.tsv`:
      ```
      chr1    1000    2000    50
      chr1    1500    2500    60
      ```
    - `output/work/*/coverage_mock_sample2.bam.tsv`:
      ```
      chr1    1000    2000    50
      chr1    1500    2500    60
      ```
  - **Why**:
    - Each BAM gets processed independently, producing identical mock coverage (in a real pipeline, coverage would vary per BAM).
    - Nextflow runs these in parallel, as seen in:
      ```
      [e9/6555a1] coverageComputation (mock_sample1.bam) | 2 of 2 ✔
      ```

#### 3. **Merge Coverage (Process: `mergeCoverage`)**
- **What It Does**:
  - Collects all coverage files (`coverage_*.tsv`) from `coverageComputation`.
  - Combines them into a single file, `merged_coverage.tsv`.
  - Simulates preparing data for querying across all BAMs.
- **Code**:
  ```groovy
  process mergeCoverage {
      input:
      path cov_files
      output:
      path "merged_coverage.tsv"
      script:
      """
      cat ${cov_files} > merged_coverage.tsv
      """
  }
  ```
  - `path cov_files`: Gets all `coverage_*.tsv` files (e.g., `coverage_mock_sample1.bam.tsv`, `coverage_mock_sample2.bam.tsv`).
  - `cat ${cov_files}`: Concatenates them into `merged_coverage.tsv`.
  - Runs after collecting outputs:
    ```groovy
    merged = mergeCoverage(coverage.collect())
    ```

- **Output**:
  - `output/work/*/merged_coverage.tsv`:
    ```
    chr1    1000    2000    50
    chr1    1500    2500    60
    chr1    1000    2000    50
    chr1    1500    2500    60
    ```
  - **Why**:
    - Combines coverage from both BAMs (four lines, two per BAM).
    - In a real pipeline, you might sum coverage or deduplicate intervals; here, it’s a simple merge.
    - Seen in:
      ```
      [a8/5660e9] mergeCoverage | 1 of 1 ✔
      ```

#### 4. **Interval Query (Process: `intervalQuery`)**
- **What It Does**:
  - Takes `merged_coverage.tsv` and one query from `queriesChannel` (e.g., `chr1:1000-2000`).
  - Runs `rust_interval_query_tool` to find coverage intervals that overlap the query.
  - Produces a TSV file with overlapping intervals for each query.
- **Code**:
  ```groovy
  process intervalQuery {
      tag "${query_interval}"
      input:
      path merged_coverage
      val(query_interval)
      output:
      path "query_result_${query_interval}.tsv"
      script:
      """
      rust_interval_query_tool \
          --interval-file "${merged_coverage}" \
          --query "${query_interval}" \
          > "query_result_${query_interval}.tsv"
      """
  }
  ```
  - `tag "${query_interval}"`: Labels tasks (e.g., `chr1:1000-2000`).
  - `path merged_coverage`: Uses `merged_coverage.tsv`.
  - `val(query_interval)`: Gets one query (e.g., `chr1:1000-2000`).
  - Runs:
    ```
    rust_interval_query_tool --interval-file "merged_coverage.tsv" --query "chr1:1000-2000" > "query_result_chr1:1000-2000.tsv"
    ```

- **Rust Tool (`rust_interval_query_tool`)**:
  - Defined in `main.rs`:
    ```rust
    #[derive(Parser, Debug)]
    struct Cli {
        #[arg(long)]
        interval_file: String,
        #[arg(long)]
        query: String,
    }
    struct Interval {
        chrom: String,
        start: u64,
        end: u64,
        coverage: u64,
    }
    fn main() -> Result<()> {
        let cli = Cli::parse();
        let file = File::open(&cli.interval_file)?;
        let reader = BufReader::new(file);
        let mut intervals = Vec::new();
        for line in reader.lines() {
            let line = line?;
            if !line.trim().is_empty() {
                intervals.push(parse_tsv_line(&line)?);
            }
        }
        let tree = IntervalTree::build(&intervals)
            .ok_or_else(|| anyhow::anyhow!("No intervals provided"))?;
        let (qchrom, qstart, qend) = parse_query(&cli.query)?;
        let hits = tree.query(qstart, qend);
        for hit in hits {
            if hit.chrom == qchrom {
                println!("{}\t{}\t{}\t{}", hit.chrom, hit.start, hit.end, hit.coverage);
            }
        }
        Ok(())
    }
    ```
  - **What It Does**:
    - Reads `--interval-file` (e.g., `merged_coverage.tsv`).
    - Parses TSV lines into `Interval` structs (e.g., `chrom: chr1, start: 1000, end: 2000, coverage: 50`).
    - Builds an interval tree for fast overlap queries.
    - Parses `--query` (e.g., `chr1:1000-2000`) into chromosome (`chr1`), start (`1000`), end (`2000`).
    - Queries the tree for intervals overlapping `qstart` to `qend` on `qchrom`.
    - Prints matching intervals in TSV format:
      ```
      chr1    1000    2000    50
      chr1    1500    2500    60
      ```

- **Output**:
  - Two files (one per query):
    - `output/work/*/query_result_chr1:1000-2000.tsv`:
      ```
      chr1    1000    2000    50
      chr1    1500    2500    60
      ```
    - `output/work/*/query_result_chr1:1500-2500.tsv`:
      ```
      chr1    1000    2000    50
      chr1    1500    2500    60
      ```
  - **Why**:
    - For `chr1:1000-2000`, it finds intervals overlapping 1000-2000 on chr1:
      - `chr1:1000-2000` (fully overlaps).
      - `chr1:1500-2500` (partially overlaps, starts at 1500).
    - Same for `chr1:1500-2500`, which overlaps both intervals.
    - Nextflow runs these in parallel:
      ```
      [5c/f4e1a8] intervalQuery (chr1:1000-2000) | 2 of 2 ✔
      ```

#### 5. **Workflow Coordination**
- **Code**:
  ```groovy
  workflow {
      coverage = coverageComputation(bamChannel)
      merged = mergeCoverage(coverage.collect())
      query_results = intervalQuery(merged, queriesChannel)
  }
  ```
- **What It Does**:
  - Runs `coverageComputation` for each BAM in `bamChannel` (2 tasks).
  - Collects all coverage files and runs `mergeCoverage` (1 task).
  - Runs `intervalQuery` for each query in `queriesChannel`, passing `merged_coverage.tsv` (2 tasks).
  - Total: 5 tasks, as seen in:
    ```
    executor >  local (5)
    ```

---

### Output Explanation
The pipeline produces files in `output/work/`, organized by task IDs (e.g., `e9/6555a1`). You can list them:
```bash
ls output/work/*/*.tsv
```

**Output Files**:
1. **`coverage_mock_sample1.bam.tsv`** and **`coverage_mock_sample2.bam.tsv`**:
   - **Content** (each):
     ```
     chr1    1000    2000    50
     chr1    1500    2500    60
     ```
   - **Meaning**:
     - Mock coverage for each BAM.
     - Shows two genomic regions on chr1 with coverage of 50 and 60 reads.
     - In a real pipeline, these would reflect actual read counts from sequencing data.

2. **`merged_coverage.tsv`**:
   - **Content**:
     ```
     chr1    1000    2000    50
     chr1    1500    2500    60
     chr1    1000    2000    50
     chr1    1500    2500    60
     ```
   - **Meaning**:
     - Combines coverage from both BAMs.
     - Has four lines (two per BAM, identical in mock data).
     - Used as input for querying to check overlaps across all samples.

3. **`query_result_chr1:1000-2000.tsv`** and **`query_result_chr1:1500-2500.tsv`**:
   - **Content** (each):
     ```
     chr1    1000    2000    50
     chr1    1500    2500    60
     ```
   - **Meaning**:
     - Lists coverage intervals from `merged_coverage.tsv` that overlap the query region.
     - For `chr1:1000-2000`:
       - `chr1:1000-2000` overlaps (starts and ends within query).
       - `chr1:1500-2500` overlaps (starts at 1500, within 1000-2000).
     - For `chr1:1500-2500`:
       - Both intervals overlap (1000-2000 starts before 2500, 1500-2500 is fully within).
     - In bioinformatics, this shows which regions have enough coverage for analysis (e.g., finding mutations).

- **Combined Output** (Optional):
  - You created `experiment_06_output.txt` manually:
    ```bash
    echo "Experiment 06 Query Results" > experiment_06_output.txt
    for file in output/work/*/query_result_*.tsv; do
        echo "Results from $file:" >> experiment_06_output.txt
        cat "$file" >> experiment_06_output.txt
        echo "" >> experiment_06_output.txt
    done
    ```
  - **Content**:
    ```
    Experiment 06 Query Results
    ---------------------------
    Results from output/work/.../query_result_chr1:1000-2000.tsv:
    chr1    1000    2000    50
    chr1    1500    2500    60

    Results from output/work/.../query_result_chr1:1500-2500.tsv:
    chr1    1000    2000    50
    chr1    1500    2500    60
    ```
  - **Meaning**:
    - Summarizes all query results in one file for easy review.
    - Useful for checking or sharing results without digging into `output/work/`.

---

### Why It Works
The corrected codes ensured success by:
- **Nextflow (`main.nf`)**:
  - Fixed DSL2 syntax (e.g., `val(bam_file)`, `path`) to avoid parsing errors.
  - Quoted arguments (e.g., `"${bam_file}"`) for robustness.
  - Trimmed channel inputs to prevent newline issues.
- **Rust Tools**:
  - `rust_coverage_tool`: Outputs mock TSV data compatible with `rust_interval_query_tool`.
  - `rust_interval_query_tool`: Reads TSV files, handles chromosomes, and outputs query results in TSV format, matching pipeline expectations.
- **Inputs**:
  - Corrected `bams.txt` to list BAM names, not commands.
  - Used valid `genome_intervals.txt` for queries.

This produced the expected files, as confirmed by your run.

---

### Real-World Context
In bioinformatics, this pipeline simulates:
- **Coverage Analysis**: Checking if sequencing data covers genomic regions well (e.g., for variant calling).
- **Interval Queries**: Finding overlaps to focus on regions of interest (e.g., genes).
- **Scalability**: Nextflow can handle thousands of BAMs and queries on a cluster, but we used mock data for learning.

The mock outputs show fake coverage, but the structure mimics tools like `samtools depth` or `bedtools coverage`.

---

### Checking Results
You can explore further:
- View a result:
  ```bash
  cat output/work/*/query_result_chr1:1000-2000.tsv
  ```
- Open in Cursor:
  ```bash
  code experiment_06_output.txt
  ```
- Count overlaps:
  ```bash
  wc -l output/work/*/query_result_*.tsv
  ```

This pipeline takes fake DNA data, calculates pretend coverage, combines it, and finds which parts overlap your regions of interest, giving you tidy TSV files to check the results.
