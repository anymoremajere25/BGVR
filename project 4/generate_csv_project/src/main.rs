use std::fs::File;
use std::io::{self, Write};
use rand::Rng;  // Used to generate random numbers for features

fn generate_csv_file() -> io::Result<()> {
    let mut file = File::create("gene_expression.csv")?;
    
    // Number of rows you want to generate
    let n_rows = 100;  // Adjust this number as needed
    let n_features = 10000;  // 10,000 features per row

    let mut rng = rand::thread_rng();

    // Generate random data for n_rows rows
    for _ in 0..n_rows {
        let features: Vec<String> = (0..n_features)
            .map(|_| rng.gen_range(0.0..10.0).to_string())  // Random floating-point values for features
            .collect();

        // Random label (either 0 or 1)
        let label: u8 = rng.gen_range(0..2);

        // Combine features and label into a single line
        let line = features.join(",") + &format!(",{}", label);

        // Write the line to the file
        writeln!(file, "{}", line)?;
    }

    Ok(())
}

fn main() -> io::Result<()> {
    generate_csv_file()?;
    println!("CSV file generated successfully!");
    Ok(())
}
