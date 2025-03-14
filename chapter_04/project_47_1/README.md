## 4.7 eQTL MAPPING AND FUNCTIONAL VARIANT DISCOVERY
### Explanation of the Rust Code for eQTL Analysis

This Rust program performs an expression quantitative trait loci (eQTL) analysis, which identifies statistical associations between genetic variants (SNPs) and gene expression levels. The program uses linear regression to model these relationships and computes p-values using the Student’s t-distribution to assess statistical significance.
### Key Components of the Code
#### 1. Data Structures

The code defines three main structs:

    Genotype: Represents SNP data across multiple individuals.

pub struct Genotype {
    pub snp_id: String,
    pub values: Vec<f64>,  // Encoded as 0, 1, or 2
}

Expression: Represents gene expression levels across the same individuals.

pub struct Expression {
    pub gene_id: String,
    pub values: Vec<f64>,
}

EqtlResult: Stores the computed slope and p-value from the linear regression for each SNP-gene pair.

    pub struct EqtlResult {
        pub snp_id: String,
        pub gene_id: String,
        pub slope: f64,
        pub p_value: f64,
    }

#### 2. Linear Regression for eQTL Analysis

The linear_eqtl function performs linear regression for a given SNP and gene expression dataset:

fn linear_eqtl(snp: &Genotype, expr: &Expression) -> EqtlResult

Steps in this function:

    Calculate means of SNP values (X̄) and expression values (Ȳ).
    Compute the slope (β1) of the linear regression:
    β1=∑(Xi−Xˉ)(Yi−Yˉ)∑(Xi−Xˉ)2
    β1​=∑(Xi​−Xˉ)2∑(Xi​−Xˉ)(Yi​−Yˉ)​
    Calculate residual variance to determine the standard error.
    Compute the Student’s t-statistic for the slope:
    t=β1SEβ1
    t=SEβ1​​β1​​ where SE_β1 = sqrt(var_resid / ss_xx).
    Compute the p-value using a Student’s t-distribution.

This tells us whether the SNP's genotype significantly affects gene expression.
#### 3. Parallelized eQTL Computation

The program processes multiple SNPs in parallel using Rayon’s parallel iterators:

let eqtl_results: Vec<EqtlResult> = all_snps
    .par_iter()
    .map(|snp| {
        genes
            .iter()
            .map(|gene| linear_eqtl(snp, gene))  // Compute for each gene
            .collect::<Vec<_>>()
    })
    .reduce_with(|mut a, mut b| { a.append(&mut b); a })
    .unwrap_or_default();

    Each SNP is processed in parallel.
    Each SNP is tested against all genes using linear_eqtl.
    The results are combined into a single list (eqtl_results).

#### 4. Writing eQTL Results to CSV

Once eQTL calculations are complete, the results are written to a CSV file:

let file = File::create("partial_eqtl.csv")?;
let mut writer = BufWriter::new(file);
for res in &eqtl_results {
    let line = format!("{},{},{:.3},{:.5}\n", res.snp_id, res.gene_id, res.slope, res.p_value);
    writer.write_all(line.as_bytes())?;
}

    Each line in the CSV file contains:

SNP_ID, GENE_ID, SLOPE, P_VALUE

Example output:

    rs1,GeneA,2.000,0.04052
    rs1,GeneB,0.150,0.43486
    rs2,GeneA,-0.750,0.41953
    rs2,GeneB,-0.125,0.34433

    The p-values indicate statistical significance.

#### 5. Running the Program

When executed, the program processes the data and prints:

Wrote 4 eQTL results.

This confirms that the analysis was completed successfully.
### Summary

    Purpose: This Rust program performs eQTL analysis using linear regression.
    Efficiency: Uses parallel processing (rayon) for faster computations.
    Statistical Methods: Computes slope, residuals, t-statistic, and p-value using Student’s t-distribution.
    Output: A CSV file containing SNP-Gene associations with statistical significance.

This approach can be extended to large genomic datasets for discovering regulatory genetic variants affecting gene expression! 
