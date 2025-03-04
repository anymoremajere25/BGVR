use reqwest::blocking;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    // Replace with a valid FASTQ file URL from NCBI FTP server
    let fastq_url = "https://ftp.ncbi.nlm.nih.gov/sra/sdk/2.10.8/samples/SRR10004850/SRR10004850.fastq";
    
    // Output path for the downloaded FASTQ
    let output_path = "downloaded_sample.fastq";

    // Create an HTTP client and issue a GET request
    let response = blocking::get(fastq_url)?;
    
    // Ensure the request was successful
    if !response.status().is_success() {
        return Err(format!("Failed to fetch FASTQ. HTTP status: {}", response.status()).into());
    }

    // Open a local file for writing
    let file = File::create(output_path)?;
    let mut writer = BufWriter::new(file);

    // Stream the response body into the local file
    let bytes = response.bytes()?;
    writer.write_all(&bytes)?;

    println!("FASTQ file has been successfully downloaded to '{}'.", output_path);

    Ok(())
}
