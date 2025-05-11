use noodles_bam::bai;
use noodles_bam::reader::Reader;
use noodles_bam::record::{Record, ReferenceSequenceId};
use noodles_core::region::Region;
use serde::Serialize;
use std::error::Error;
use std::fs::File;
use std::io::{self, BufReader};
use std::path::Path;
use std::ops::Bound;

#[derive(Serialize, Debug)]
struct RegionCoverage {
    reference_name: String,
    start: i32,
    end: i32,
    coverage: Vec<u32>,
}

fn main() -> Result<(), Box<dyn Error>> {
    let bam_path = "test.bam";
    let bai_path = format!("{}.bai", bam_path);
    let reference_name = "chr1";
    let region_start = 10000;
    let region_end = 10100;
    
    // Create region using the mapped constructor for noodles-core 0.5
    let region = Region::mapped(reference_name, region_start..=region_end);

    if !Path::new(&bai_path).exists() {
        // We'll defer this to manual indexing or use a bai file that exists
        return Err("BAI index file does not exist. Please generate it first.".into());
    }

    let coverage = compute_coverage(bam_path, &bai_path, &region)?;
    println!("{}", serde_json::to_string_pretty(&coverage)?);

    Ok(())
}

fn compute_coverage(
    bam_path: &str,
    bai_path: &str,
    region: &Region,
) -> Result<RegionCoverage, Box<dyn Error>> {
    // Load the BAM file
    let mut reader = Reader::new(BufReader::new(File::open(bam_path)?));
    
    // Read header and reference sequences (discard header as it's not needed)
    reader.read_header()?;
    let reference_sequences = reader.read_reference_sequences()?;
    
    // Read the BAI index (just check it exists - we're not using it for random access in this example)
    // We could remove this line if we're not using the index, but keeping for validation
    bai::read(bai_path)?;
    
    // Get reference sequence ID and name
    let reference_name = region.name();
    let ref_id = reference_sequences
        .get_index_of(reference_name)
        .ok_or_else(|| io::Error::new(io::ErrorKind::NotFound, "Reference not found"))?;
    
    // Create a ReferenceSequenceId using the correct method
    let ref_seq_id = ReferenceSequenceId::try_from(ref_id as i32)
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Invalid reference sequence ID"))?;
    
    // Extract start and end from the region using noodles-core 0.5 API
    let mapped = match region {
        Region::Mapped(mapped) => mapped,
        _ => return Err("Only mapped regions are supported".into()),
    };
    
    // Extract bounds directly from the interval tuple
    let interval = mapped.interval();
    
    // Direct destructuring of the tuple
    let (start_bound, end_bound) = interval;
    
    let start = match start_bound {
        Bound::Included(pos) => pos,
        Bound::Excluded(pos) => pos + 1,
        Bound::Unbounded => 0,
    };
    
    let end = match end_bound {
        Bound::Included(pos) => pos,
        Bound::Excluded(pos) => pos - 1,
        Bound::Unbounded => reference_sequences[ref_id].len() as i32,
    };
    
    // Initialize coverage array
    let region_length = (end - start + 1) as usize;
    let mut coverage = vec![0u32; region_length];
    
    // Process each record in the BAM file
    let mut record = Record::default();
    
    // Reset the reader position
    let mut reader = Reader::new(BufReader::new(File::open(bam_path)?));
    reader.read_header()?;
    reader.read_reference_sequences()?;
    
    while reader.read_record(&mut record)? != 0 {
        // Check if record matches our region by reference ID
        if let Some(record_ref_id) = record.reference_sequence_id() {
            // Now compare the correct types
            if record_ref_id == ref_seq_id {
                // Check position
                if let Some(position) = record.position() {
                    let record_start = i32::from(position) - 1; // Convert to 0-based
                    
                    // Determine end position using CIGAR
                    let cigar = record.cigar();
                    
                    // Calculate alignment length manually from CIGAR string
                    let mut alignment_length = 0;
                    
                    // Get the CIGAR string representation and parse it
                    let cigar_str = cigar.to_string();
                    let mut num_buffer = String::new();
                    
                    for c in cigar_str.chars() {
                        if c.is_ascii_digit() {
                            num_buffer.push(c);
                        } else {
                            // Found an operation character
                            let length = num_buffer.parse::<i32>().unwrap_or(0);
                            num_buffer.clear();
                            
                            // Count operations that consume reference sequence
                            match c {
                                'M' | '=' | 'X' | 'D' | 'N' => alignment_length += length,
                                _ => {}
                            }
                        }
                    }
                    
                    let record_end = record_start + alignment_length - 1;
                    
                    // Check for overlap with our region
                    if !(record_end < start || record_start > end) {
                        // Calculate overlap with target region
                        let overlap_start = record_start.max(start);
                        let overlap_end = record_end.min(end);
                        
                        // Update coverage for each position
                        for pos in overlap_start..=overlap_end {
                            let index = (pos - start) as usize;
                            if index < coverage.len() {
                                coverage[index] += 1;
                            }
                        }
                    }
                }
            }
        }
    }
    
    Ok(RegionCoverage {
        reference_name: reference_name.to_string(),
        start,
        end,
        coverage,
    })
}