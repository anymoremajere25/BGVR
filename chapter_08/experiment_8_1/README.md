**8.1. Fundamentals of Genetic Variation**
**experiment_8_1**

The following Rust code demonstrates how to compute genotype frequencies, perform a chi-square–based Hardy-Weinberg (HW) p-value analysis, and leverage concurrency to process genomic data in parallel. This example has been optimized for large-scale applications by integrating crates that support stable numerical computations, efficient data processing, and multithreading.

In this implementation:

* **rust-htslib** is used to parse VCF/BCF files.
* **rayon** enables parallel processing of genomic data across chunks.
* **ndarray** handles numerical operations needed for genotype frequency calculations.
* **statrs** provides statistical functions, including chi-square p-value computations.
* **polars** manages data as DataFrames and allows fast in-memory processing.

By adjusting memory limits, chunk sizes, and concurrency settings, this setup scales efficiently to handle large datasets typical of pharmaceutical pipelines or consortium-scale research.

**Project structure:**

```
experiment_8_1/
├── Cargo.toml                  # Rust dependencies
└── src/
    ├── main.rs                # Main Rust program
    ├── synthetic.vcf          # Input VCF file
    ├── synthetic.vcf.hw_results.csv  # Output CSV with HW p-values
    └── output.txt             # Additional output
```

**To run the program using WSL:**

```bash
cargo run -- synthetic.vcf 0 1000000
```

This will generate the file `synthetic.vcf.hw_results.csv`.

**Dependencies:**

```toml
[dependencies]
rust-htslib = "0.49.0"
rayon = "1.5.1"
ndarray = "0.16.1"
statrs = "0.18.0"
polars = { version = "0.46", features = ["lazy"] }
```

---

### 🧠 What the Code Does (main.rs)

This program reads a VCF file, analyzes each variant's genotype distribution, calculates a p-value using the Hardy-Weinberg chi-square test, and outputs the results to a CSV.

---

### 🔍 How It Works

1. **Libraries Used:**

   * `statrs`: Performs chi-square test and calculates p-values.
   * `polars`: Manages data in tabular format and handles CSV export.
   * `rust-htslib`: Mentioned but not used directly (VCF parsing is manual).

2. **Main Functions:**

**(A) `chi_square_hw(aa, ab, bb, p)`**

* Inputs:

  * `aa`, `ab`, `bb`: Counts of 0/0, 0/1, and 1/1 genotypes.
  * `p`: Reference allele frequency.
* Computes expected genotype counts under HW equilibrium.
* Calculates the chi-square statistic and corresponding p-value.

  * High p-value ≈ fits HW equilibrium.
  * Low p-value ≈ deviation from HW.

**(B) `process_vcf_file(vcf_path, start_pos, end_pos)`**

* Reads the VCF file line-by-line.
* Skips metadata (`##` lines).
* Extracts sample genotype columns from the `#CHROM` header.
* For each variant:

  * Filters by position.
  * Extracts and counts genotypes.
  * Computes allele frequency and HW p-value.
  * Stores data in a list, which is then turned into a DataFrame.
* Returns the DataFrame.

**(C) `main()`**

* Reads command-line arguments: VCF file, start and end positions.
* Calls `process_vcf_file`.
* Prints and writes results to `*.hw_results.csv`.

---

### 📂 Sample Input: `synthetic.vcf`

```vcf
##fileformat=VCFv4.2
##contig=<ID=1,length=249250621>
#CHROM  POS     ID   REF   ALT   QUAL  FILTER INFO FORMAT Sample1 Sample2 Sample3
1       12345   .    A     G     50.0  PASS   NS=3  GT     0/0     0/1     1/1
1       67890   .    T     C     40.0  PASS   NS=3  GT     0/1     1/1     0/0
```

**Details for Variant at Position 12345:**

* Genotypes: 0/0, 0/1, 1/1 → Counts: 1, 1, 1
* Allele Frequency `p` = (2×1 + 1)/(2×3) = 3/6 = 0.5

**Expected Genotype Frequencies (under HW):**

* AA: (0.5²) × 3 = 0.75
* AB: 2×0.5×0.5 × 3 = 1.5
* BB: (0.5²) × 3 = 0.75

**Chi-square Statistic:**

```text
χ² = ((1−0.75)²/0.75) + ((1−1.5)²/1.5) + ((1−0.75)²/0.75)
    = 0.3333
```

**P-value:**

```text
p = 1 - CDF(χ² = 0.3333, df = 1) ≈ 0.5637
```

**Interpretation:**

* p-value > 0.05 → No significant deviation from HW equilibrium.

---

### 📝 Output: `synthetic.vcf.hw_results.csv`

| Chromosome | Position | Reference Allele | Alternate Allele | HWE p-value |
| ---------- | -------- | ---------------- | ---------------- | ----------- |
| 1          | 12345    | A                | G                | 0.5637      |
| 1          | 67890    | T                | C                | 0.5637      |

---

### 🛠 Executable: `vcf_analysis`

Once compiled, the binary can be run like this:

```bash
./vcf_analysis synthetic.vcf
```

It will automatically generate and display the HW p-values in a CSV.

---

### ✅ Summary

* **Goal**: Check if genetic variants conform to Hardy-Weinberg equilibrium.
* **Result**: Both variants have p-values ≈ 0.56 → no significant deviation from HW.
* **Why it matters**: Deviations from HW can reveal issues such as genotyping errors, selection, or population structure.
