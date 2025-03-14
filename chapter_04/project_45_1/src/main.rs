use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};

/// Represents an exon or a segment in the splicing graph.
///
/// The fields are prefixed with `_` to avoid compiler warnings
/// if they are currently unused.
#[derive(Debug, Clone)]
struct ExonSegment {
    _start: usize,
    _end: usize,
}

/// Represents alignment data for one read.
///
/// The field is prefixed with `_` to avoid compiler warnings
/// if it is currently unused.
#[derive(Debug, Clone)]
struct Alignment {
    _chrom: String,
    start: usize,
    end: usize,
}

/// Graph-like structure holding exons and splice junctions.
/// - adjacency: a map from exon_key (start, end) to a vector of (ExonSegment, coverage)
#[derive(Debug, Default, Clone)]
struct SplicingGraph {
    adjacency: HashMap<(usize, usize), Vec<(ExonSegment, u64)>>,
}

impl SplicingGraph {
    /// Adds or updates a junction from an exon (exon_key) to a target exon.
    fn add_junction(&mut self, exon_key: (usize, usize), target_exon: ExonSegment) {
        let coverage = 1;
        self.adjacency
            .entry(exon_key)
            .or_default()
            .push((target_exon, coverage));
    }

    /// Merges another SplicingGraph into the current one.
    fn merge(&mut self, other: SplicingGraph) {
        for (key, edges) in other.adjacency {
            self.adjacency.entry(key).or_default().extend(edges);
        }
    }
}

/// Processes a slice of `Alignment` entries to produce a local splicing graph.
fn process_alignment_chunk(batch: &[Alignment]) -> SplicingGraph {
    let mut local_graph = SplicingGraph::default();
    for align in batch {
        // We do not use `_chrom`, but we do use (align.start, align.end) as exon_key
        let exon_key = (align.start, align.end);
        let target_exon = ExonSegment {
            _start: align.start + 50,
            _end: align.end + 100,
        };
        local_graph.add_junction(exon_key, target_exon);
    }
    local_graph
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Hypothetical alignment data
    let all_alignments = vec![
        Alignment {
            _chrom: "chr1".to_string(),
            start: 100,
            end: 200,
        },
        Alignment {
            _chrom: "chr1".to_string(),
            start: 150,
            end: 300,
        },
        Alignment {
            _chrom: "chr2".to_string(),
            start: 500,
            end: 700,
        },
        // ...
    ];

    // Define a chunk size
    let chunk_size = 100;

    // Partition the alignments into chunks
    let chunks: Vec<_> = all_alignments
        .chunks(chunk_size)
        .map(|c| c.to_vec())
        .collect();

    // Build and merge partial splicing graphs in parallel
    let final_graph = chunks
        .into_par_iter()
        .map(|batch| process_alignment_chunk(&batch))
        .reduce(
            SplicingGraph::default,
            |mut acc, local_graph| {
                acc.merge(local_graph);
                acc
            },
        );

    // Serialize or write the final splicing graph to disk
    let out_file = File::create("partial_splicing_graph.bin")?;
    let mut writer = BufWriter::new(out_file);
    writer.write_all(b"Serialized splicing graph example\n")?;
    writer.write_all(format!("{:#?}", final_graph).as_bytes())?;

    println!("Splicing graph has been successfully written.");

    Ok(())
}
