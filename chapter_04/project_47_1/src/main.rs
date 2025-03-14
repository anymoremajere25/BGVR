use rayon::prelude::*;
use statrs::distribution::{StudentsT, ContinuousCDF};
use std::fs::File;
use std::io::{BufWriter, Write};

/// Represents genotype data for one SNP across multiple individuals.
#[derive(Clone)]
pub struct Genotype {
    pub snp_id: String,
    pub values: Vec<f64>,  // e.g. coded 0,1,2
}

/// Expression data for one gene across the same individuals.
#[derive(Clone)]
pub struct Expression {
    pub gene_id: String,
    pub values: Vec<f64>,
}

/// Holds the result of a single eQTL test.
#[derive(Debug)]
pub struct EqtlResult {
    pub snp_id: String,
    pub gene_id: String,
    pub slope: f64,
    pub p_value: f64,
}

/// A basic linear regression: Y ~ X + error.
/// Uses a Student’s t-distribution for p-values, with n - 2 degrees of freedom.
fn linear_eqtl(snp: &Genotype, expr: &Expression) -> EqtlResult {
    assert_eq!(snp.values.len(), expr.values.len());
    let x = &snp.values;
    let y = &expr.values;
    let n = x.len() as f64;
    let df = n - 2.0; // degrees of freedom for slope estimate

    // Compute means
    let mean_x = x.iter().sum::<f64>() / n;
    let mean_y = y.iter().sum::<f64>() / n;

    // Compute slope
    let mut ss_xy = 0.0;
    let mut ss_xx = 0.0;
    for i in 0..x.len() {
        let dx = x[i] - mean_x;
        let dy = y[i] - mean_y;
        ss_xy += dx * dy;
        ss_xx += dx * dx;
    }
    let slope = ss_xy / ss_xx;

    // Compute variance of residuals
    let mut rss = 0.0; // residual sum of squares
    for i in 0..x.len() {
        let pred = mean_y + slope * (x[i] - mean_x);
        let resid = y[i] - pred;
        rss += resid * resid;
    }
    let var_resid = rss / df;             // error variance
    let se_slope = (var_resid / ss_xx).sqrt(); // std error of slope

    // Student's t-statistic for slope
    let t_stat = slope / se_slope;
    let t_dist = StudentsT::new(0.0, 1.0, df).unwrap();
    let p_value = 2.0 * (1.0 - t_dist.cdf(t_stat.abs()));

    EqtlResult {
        snp_id: snp.snp_id.clone(),
        gene_id: expr.gene_id.clone(),
        slope,
        p_value,
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Hypothetical SNP and gene data (small example).
    // In real usage, you may have millions of SNPs and thousands of genes.
    let all_snps = vec![
        Genotype {
            snp_id: "rs1".to_string(),
            values: vec![0.0, 1.0, 2.0, 1.0, 1.0],
        },
        Genotype {
            snp_id: "rs2".to_string(),
            values: vec![2.0, 2.0, 0.0, 0.0, 1.0],
        },
    ];

    let genes = vec![
        Expression {
            gene_id: "GeneA".to_string(),
            values: vec![10.0, 12.0, 14.0, 11.0, 13.0],
        },
        Expression {
            gene_id: "GeneB".to_string(),
            values: vec![8.2, 7.9, 8.5, 8.1, 8.0],
        },
    ];

    // We process each SNP in parallel, computing eQTL results for all genes locally,
    // then merge those local vectors of results into a single final vector.
    let eqtl_results: Vec<EqtlResult> = all_snps
        .par_iter()
        .map(|snp| {
            genes
                .iter()
                .map(|gene| linear_eqtl(snp, gene))
                .collect::<Vec<_>>()
        })
        .reduce_with(|mut a, mut b| {
            a.append(&mut b);
            a
        })
        .unwrap_or_default();

    // Write results to a CSV file
    let file = File::create("partial_eqtl.csv")?;
    let mut writer = BufWriter::new(file);
    for res in &eqtl_results {
        let line = format!("{},{},{:.3},{:.5}\n", res.snp_id, res.gene_id, res.slope, res.p_value);
        writer.write_all(line.as_bytes())?;
    }

    println!("Wrote {} eQTL results.", eqtl_results.len());
    Ok(())
}
