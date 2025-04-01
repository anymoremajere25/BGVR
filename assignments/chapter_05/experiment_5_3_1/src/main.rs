use anyhow::{Context, Result};
use rayon::prelude::*;
use bio::data_structures::suffix_array::suffix_array;
use bio::data_structures::bwt::{bwt, Less, Occ};
use bio::data_structures::fmindex::{FMIndex, FMIndexable};
use bio::alphabets::Alphabet;
use serde::{Deserialize, Serialize};
use std::fs::{File};
use std::io::{BufWriter};
use std::path::PathBuf;
use clap::Parser;

/// Parameters for FM-index-based alignment.
#[derive(Parser, Debug)]
struct Args {
    #[arg(long)]
    reference: PathBuf,

    #[arg(long)]
    reads: PathBuf,

    #[arg(long, default_value = "fmindex_chunks")]
    index_outdir: PathBuf,

    #[arg(long, default_value_t = 50)]
    reference_chunk_size: usize,

    #[arg(long, default_value_t = 3)]
    fm_sampling_rate: usize,

    #[arg(long, default_value = "alignments.json")]
    alignment_output: PathBuf,
}

#[derive(Serialize, Deserialize, Debug)]
struct ReadAlignment {
    read_seq: Vec<u8>,
    matches: Vec<usize>,
}

/// Preprocess reference genome: ensure uppercase & remove invalid characters
fn preprocess_reference(reference: &[u8]) -> Vec<u8> {
    let mut processed: Vec<u8> = reference
        .iter()
        .filter(|&&c| b"ACGTN".contains(&c)) // Keep only valid nucleotides
        .map(|&c| c.to_ascii_uppercase()) // Ensure uppercase
        .collect();
    
    if !processed.ends_with(&[b'$']) {
        processed.push(b'$'); // Append sentinel character if missing
    }

    processed
}

/// Build the FM-index from the reference genome
fn build_fmindex(reference: &[u8], sampling_rate: usize) -> Result<(FMIndex<Vec<u8>, Vec<usize>, Occ>, Vec<usize>)> {
    if reference.is_empty() {
        return Err(anyhow::anyhow!("❌ Reference genome is empty!"));
    }

    println!("🔹 Reference after preprocessing: {:?}", String::from_utf8_lossy(&reference));

    let sa = suffix_array(reference);
    if sa.is_empty() {
        return Err(anyhow::anyhow!("❌ Suffix array is empty!"));
    }

    println!("✅ Suffix array length: {}", sa.len());
    println!("🔹 First 10 suffix array indices: {:?}", &sa[..10.min(sa.len())]);

    let bwt_data = bwt(reference, &sa);
    if bwt_data.is_empty() {
        return Err(anyhow::anyhow!("❌ BWT computation failed! BWT is empty."));
    }

    println!("✅ BWT length: {}", bwt_data.len());
    println!("🔹 First 10 BWT characters: {:?}", String::from_utf8_lossy(&bwt_data[..10.min(bwt_data.len())]));

    let alphabet = Alphabet::new(b"ACGTN");
    let less = Less::new();
    let occ = Occ::new(&bwt_data, sampling_rate as u32, &alphabet);

    Ok((FMIndex::new(bwt_data, less, occ), sa))
}

/// Perform parallel alignment of reads against the FM-index
fn parallel_align_reads(
    reference: &[u8],
    _reference_chunk_size: usize,
    reads: &[Vec<u8>],
) -> Result<Vec<ReadAlignment>> {
    let (fm, sa) = build_fmindex(reference, 3)?;

    let alignments: Vec<ReadAlignment> = reads
        .par_iter()
        .map(|read| {
            println!("🔍 Searching for read: {:?}", String::from_utf8_lossy(read));

            match fm.backward_search(read.iter()) {
                bio::data_structures::fmindex::BackwardSearchResult::Complete(interval) => {
                    println!("✅ Interval found: {:?}", interval);

                    if interval.upper > sa.len() {
                        eprintln!("⚠ Warning: interval.upper ({}) exceeds suffix array length ({})", interval.upper, sa.len());
                        return ReadAlignment {
                            read_seq: read.clone(),
                            matches: vec![],
                        };
                    }

                    if interval.lower < interval.upper {
                        let positions: Vec<usize> = sa[interval.lower..interval.upper].to_vec();
                        println!("✅ Found {} matches for read", positions.len());

                        ReadAlignment {
                            read_seq: read.clone(),
                            matches: positions,
                        }
                    } else {
                        eprintln!("⚠ Warning: Invalid interval found ({:?})", interval);
                        ReadAlignment {
                            read_seq: read.clone(),
                            matches: vec![],
                        }
                    }
                }
                _ => {
                    println!("❌ No match found for read.");
                    ReadAlignment {
                        read_seq: read.clone(),
                        matches: vec![],
                    }
                }
            }
        })
        .collect();

    Ok(alignments)
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Read and preprocess reference genome
    let reference_raw = std::fs::read(&args.reference)
        .with_context(|| format!("❌ Failed to read reference genome from {:?}", args.reference))?;
    let reference = preprocess_reference(&reference_raw);

    println!("✅ Processed reference genome length: {}", reference.len());

    // Read reads file
    let reads_data = std::fs::read_to_string(&args.reads)
        .with_context(|| format!("❌ Failed to read reads file from {:?}", args.reads))?;
    let reads: Vec<Vec<u8>> = reads_data.lines().map(|line| line.as_bytes().to_vec()).collect();

    println!("✅ Loaded {} reads", reads.len());

    // Perform read alignment
    let alignments = parallel_align_reads(&reference, args.reference_chunk_size, &reads)?;

    // Write alignments to file
    let file = File::create(&args.alignment_output)
        .with_context(|| format!("❌ Failed to create alignment output file {:?}", args.alignment_output))?;
    serde_json::to_writer(&mut BufWriter::new(file), &alignments)
        .with_context(|| format!("❌ Failed to write alignment results to {:?}", args.alignment_output))?;

    println!("✅ Wrote alignment results for {} reads to {:?}", alignments.len(), args.alignment_output);

    Ok(())
}