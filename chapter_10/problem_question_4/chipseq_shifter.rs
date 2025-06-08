use std::cmp;

#[derive(Debug, Clone, PartialEq)]
pub struct Read {
    pub chromosome: String,
    pub start: u64,
    pub end: u64,
    pub strand: char, // '+' for forward, '-' for reverse
}

impl Read {
    pub fn new(chromosome: String, start: u64, end: u64, strand: char) -> Self {
        Read {
            chromosome,
            start,
            end,
            strand,
        }
    }
}

/// Shifts ChIP-Seq reads to better approximate transcription factor binding sites.
/// 
/// # Arguments
/// * `reads` - Vector of Read structs representing aligned sequencing reads
/// * `shift_size` - Optional shift distance in base pairs (default: 100)
/// 
/// # Returns
/// Vector of Read structs with adjusted coordinates
/// 
/// # Behavior
/// - Forward strand reads (+): shifted downstream (coordinates increased)
/// - Reverse strand reads (-): shifted upstream (coordinates decreased)
/// - Handles edge cases like negative coordinates and potential overflow
pub fn shift_reads(reads: Vec<Read>, shift_size: Option<u32>) -> Vec<Read> {
    let shift = shift_size.unwrap_or(100) as u64;
    
    reads.into_iter().map(|read| {
        shift_single_read(read, shift)
    }).collect()
}

/// Shifts a single read based on its strand orientation
fn shift_single_read(mut read: Read, shift: u64) -> Read {
    let read_length = read.end - read.start;
    
    match read.strand {
        '+' => {
            // Forward strand: shift downstream (increase coordinates)
            // Check for potential overflow
            if read.start <= u64::MAX - shift {
                read.start += shift;
                read.end = read.start + read_length;
            } else {
                // Handle overflow by capping at maximum value
                read.start = u64::MAX - read_length;
                read.end = u64::MAX;
            }
        },
        '-' => {
            // Reverse strand: shift upstream (decrease coordinates)
            if read.start >= shift {
                read.start -= shift;
                read.end = read.start + read_length;
            } else {
                // Handle underflow by setting to start of chromosome
                read.start = 0;
                read.end = read_length;
            }
        },
        _ => {
            // Invalid strand character - return original read unchanged
            eprintln!("Warning: Invalid strand '{}' for read at {}:{}-{}", 
                     read.strand, read.chromosome, read.start, read.end);
        }
    }
    
    read
}

/// Concurrent version for large datasets using rayon
#[cfg(feature = "parallel")]
pub fn shift_reads_parallel(reads: Vec<Read>, shift_size: Option<u32>) -> Vec<Read> {
    use rayon::prelude::*;
    
    let shift = shift_size.unwrap_or(100) as u64;
    
    reads.into_par_iter().map(|read| {
        shift_single_read(read, shift)
    }).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_forward_strand_shift() {
        let reads = vec![
            Read::new("chr1".to_string(), 1000, 1050, '+'),
            Read::new("chr2".to_string(), 2000, 2075, '+'),
        ];
        
        let shifted = shift_reads(reads, Some(100));
        
        assert_eq!(shifted[0].start, 1100);
        assert_eq!(shifted[0].end, 1150);
        assert_eq!(shifted[1].start, 2100);
        assert_eq!(shifted[1].end, 2175);
    }

    #[test]
    fn test_reverse_strand_shift() {
        let reads = vec![
            Read::new("chr1".to_string(), 1000, 1050, '-'),
            Read::new("chr2".to_string(), 2000, 2075, '-'),
        ];
        
        let shifted = shift_reads(reads, Some(100));
        
        assert_eq!(shifted[0].start, 900);
        assert_eq!(shifted[0].end, 950);
        assert_eq!(shifted[1].start, 1900);
        assert_eq!(shifted[1].end, 1975);
    }

    #[test]
    fn test_default_shift_size() {
        let reads = vec![
            Read::new("chr1".to_string(), 1000, 1050, '+'),
        ];
        
        let shifted = shift_reads(reads, None);
        
        assert_eq!(shifted[0].start, 1100); // 1000 + 100 (default)
        assert_eq!(shifted[0].end, 1150);
    }

    #[test]
    fn test_underflow_protection() {
        let reads = vec![
            Read::new("chr1".to_string(), 50, 100, '-'), // Would underflow with 100bp shift
        ];
        
        let shifted = shift_reads(reads, Some(100));
        
        assert_eq!(shifted[0].start, 0);
        assert_eq!(shifted[0].end, 50);
    }

    #[test]
    fn test_overflow_protection() {
        let reads = vec![
            Read::new("chr1".to_string(), u64::MAX - 25, u64::MAX, '+'),
        ];
        
        let shifted = shift_reads(reads, Some(100));
        
        // Should be capped to prevent overflow
        assert_eq!(shifted[0].start, u64::MAX - 25);
        assert_eq!(shifted[0].end, u64::MAX);
    }

    #[test]
    fn test_invalid_strand() {
        let reads = vec![
            Read::new("chr1".to_string(), 1000, 1050, 'x'), // Invalid strand
        ];
        
        let shifted = shift_reads(reads, Some(100));
        
        // Should remain unchanged
        assert_eq!(shifted[0].start, 1000);
        assert_eq!(shifted[0].end, 1050);
    }

    #[test]
    fn test_mixed_strands() {
        let reads = vec![
            Read::new("chr1".to_string(), 1000, 1050, '+'),
            Read::new("chr1".to_string(), 2000, 2050, '-'),
            Read::new("chr2".to_string(), 3000, 3100, '+'),
        ];
        
        let shifted = shift_reads(reads, Some(75));
        
        // Forward strand: +75
        assert_eq!(shifted[0].start, 1075);
        assert_eq!(shifted[0].end, 1125);
        
        // Reverse strand: -75
        assert_eq!(shifted[1].start, 1925);
        assert_eq!(shifted[1].end, 1975);
        
        // Forward strand: +75
        assert_eq!(shifted[2].start, 3075);
        assert_eq!(shifted[2].end, 3175);
    }

    #[test]
    fn test_preserve_read_length() {
        let reads = vec![
            Read::new("chr1".to_string(), 1000, 1150, '+'), // 150bp read
            Read::new("chr1".to_string(), 2000, 2025, '-'), // 25bp read
        ];
        
        let shifted = shift_reads(reads.clone(), Some(50));
        
        // Read lengths should be preserved
        assert_eq!(shifted[0].end - shifted[0].start, reads[0].end - reads[0].start);
        assert_eq!(shifted[1].end - shifted[1].start, reads[1].end - reads[1].start);
    }
}

// Example usage and demonstration
fn main() {
    println!("ChIP-Seq Read Shifter Demo");
    println!("==========================");
    
    // Create synthetic ChIP-Seq reads
    let reads = vec![
        Read::new("chr1".to_string(), 1000, 1050, '+'),
        Read::new("chr1".to_string(), 1200, 1250, '-'),
        Read::new("chr2".to_string(), 5000, 5075, '+'),
        Read::new("chr2".to_string(), 5500, 5550, '-'),
        Read::new("chrX".to_string(), 100, 150, '+'),
        Read::new("chrX".to_string(), 200, 250, '-'),
    ];
    
    println!("Original reads:");
    for (i, read) in reads.iter().enumerate() {
        println!("  Read {}: {}:{}-{} ({})", 
                i + 1, read.chromosome, read.start, read.end, read.strand);
    }
    
    // Shift reads with default 100bp
    let shifted_default = shift_reads(reads.clone(), None);
    println!("\nAfter shifting by 100bp (default):");
    for (i, read) in shifted_default.iter().enumerate() {
        println!("  Read {}: {}:{}-{} ({})", 
                i + 1, read.chromosome, read.start, read.end, read.strand);
    }
    
    // Shift reads with custom 150bp
    let shifted_custom = shift_reads(reads.clone(), Some(150));
    println!("\nAfter shifting by 150bp:");
    for (i, read) in shifted_custom.iter().enumerate() {
        println!("  Read {}: {}:{}-{} ({})", 
                i + 1, read.chromosome, read.start, read.end, read.strand);
    }
    
    // Demonstrate edge case handling
    let edge_case_reads = vec![
        Read::new("chr1".to_string(), 50, 100, '-'),  // Will underflow
        Read::new("chr1".to_string(), u64::MAX - 25, u64::MAX, '+'), // Will overflow
    ];
    
    println!("\nEdge case handling:");
    println!("Before: chr1:50-100 (-), chr1:{}-{} (+)", 
             u64::MAX - 25, u64::MAX);
    
    let shifted_edge = shift_reads(edge_case_reads, Some(100));
    println!("After:  chr1:{}-{} (-), chr1:{}-{} (+)", 
             shifted_edge[0].start, shifted_edge[0].end,
             shifted_edge[1].start, shifted_edge[1].end);
}