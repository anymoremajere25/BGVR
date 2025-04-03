use anyhow::{Context, Result};
use rayon::prelude::*;
use needletail::{parse_fastx_file, FastxReader};
use fnv::FnvHashMap;
use serde::{Serialize, Deserialize};
use std::{
    fs::{File, create_dir_all},
    io::{BufWriter, BufReader, Write},
    path::PathBuf,
};
use clap::Parser;

/// Command-line arguments for chunk-based k-mer counting and de Bruijn graph construction.
#[derive(Parser, Debug)]
#[command(name = "kmer_debruijn_builder")]
#[command(about = "Performs chunked k-mer counting and builds a minimal de Bruijn graph")]
struct Args {
    /// The input FASTQ file path (supports gzipped FASTQ).
    #[arg(long)]
    input: PathBuf,

    /// The size of k.
    #[arg(long, default_value_t = 31)]
    k: usize,

    /// Minimum count threshold to include a k-mer in the de Bruijn graph.
    #[arg(long, default_value_t = 2)]
    threshold: u64,

    /// Number of sequences to read per chunk, to avoid loading the entire FASTQ at once.
    #[arg(long, default_value_t = 10_000)]
    chunk_size: usize,

    /// Directory in which to store partial k-mer maps from each chunk.
    #[arg(long, default_value = "partial_kmer_maps")]
    partial_outdir: PathBuf,

    /// Name of the final merged adjacency file.
    #[arg(long, default_value = "final_debruijn.bin")]
    final_output: PathBuf,
}

/// A minimal de Bruijn graph structure. Each entry in `adjacency` tracks
/// a prefix (k-1-mer) and a map of next base -> count.
#[derive(Serialize, Deserialize, Debug, Default)]
struct DeBruijn {
    adjacency: FnvHashMap<Vec<u8>, FnvHashMap<u8, u64>>,
}

impl DeBruijn {
    /// Creates an empty de Bruijn graph.
    fn new() -> Self {
        Self {
            adjacency: FnvHashMap::default(),
        }
    }

    /// Incorporates a full k-mer into the graph by splitting it into a (k-1)-mer prefix and next base.
    fn add_kmer(&mut self, kmer: &[u8]) {
        if kmer.len() > 1 {
            let prefix = &kmer[..kmer.len() - 1];
            let next_base = kmer[kmer.len() - 1];
            let entry = self.adjacency.entry(prefix.to_vec()).or_default();
            *entry.entry(next_base).or_insert(0) += 1;
        }
    }

    /// Merges another de Bruijn graph into this one by aggregating edge counts.
    fn merge(&mut self, other: DeBruijn) {
        for (prefix, edges) in other.adjacency {
            let main_entry = self.adjacency.entry(prefix).or_default();
            for (base, count) in edges {
                *main_entry.entry(base).or_insert(0) += count;
            }
        }
    }
}

/// A function to build a de Bruijn graph from a k-mer map above a certain count threshold.
fn build_debruijn(kmer_map: &FnvHashMap<Vec<u8>, u64>, threshold: u64) -> DeBruijn {
    let mut dbg = DeBruijn::new();
    for (kmer, &count) in kmer_map.iter() {
        if count >= threshold {
            dbg.add_kmer(kmer);
        }
    }
    dbg
}

/// Serializes a de Bruijn graph to disk in a simple binary format (via bincode).
fn write_debruijn_graph(dbg: &DeBruijn, path: &PathBuf) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("Failed to create de Bruijn output file at {:?}", path))?;
    bincode::serialize_into(BufWriter::new(file), dbg)
        .with_context(|| format!("Failed to serialize de Bruijn graph to {:?}", path))?;
    Ok(())
}

/// Deserializes a de Bruijn graph from disk.
fn read_debruijn_graph(path: &PathBuf) -> Result<DeBruijn> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open de Bruijn file at {:?}", path))?;
    let dbg: DeBruijn = bincode::deserialize_from(BufReader::new(file))
        .with_context(|| format!("Failed to deserialize de Bruijn graph from {:?}", path))?;
    Ok(dbg)
}

/// K-mer counting for a slice of sequences, returning a local map from k-mer -> count.
fn count_kmers_in_records(records: &[needletail::FastxRecord], k: usize) -> FnvHashMap<Vec<u8>, u64> {
    let mut local_map = FnvHashMap::default();
    for seqrec in records {
        let seq = seqrec.seq();
        for i in 0..=seq.len().saturating_sub(k) {
            let kmer = &seq[i..(i + k)];
            *local_map.entry(kmer.to_vec()).or_insert(0) += 1;
        }
    }
    local_map
}

/// Merges two k-mer maps.
fn merge_kmer_maps(mut global_map: FnvHashMap<Vec<u8>, u64>, local_map: FnvHashMap<Vec<u8>, u64>) -> FnvHashMap<Vec<u8>, u64> {
    for (kmer, count) in local_map {
        *global_map.entry(kmer).or_insert(0) += count;
    }
    global_map
}

/// Reads a chunk of records from the FASTQ reader, up to chunk_size.
fn read_chunk<R: FastxReader>(reader: &mut R, chunk_size: usize) -> Result<Vec<needletail::FastxRecord>> {
    let mut chunk = Vec::with_capacity(chunk_size);
    for _ in 0..chunk_size {
        if let Some(res) = reader.next() {
            let record = res?;
            chunk.push(record);
        } else {
            break;
        }
    }
    Ok(chunk)
}

fn main() -> Result<()> {
    let args = Args::parse();

    create_dir_all(&args.partial_outdir)
        .with_context(|| format!("Failed to create partial output directory at {:?}", args.partial_outdir))?;

    // Open the FASTQ file.
    let mut reader = parse_fastx_file(&args.input)
        .with_context(|| format!("Failed to open FASTQ file at {:?}", args.input))?;

    // In a real HPC environment, ephemeral tasks could each handle one or more chunks.
    // Here, we demonstrate a single process reading chunks sequentially, building partial k-mer maps.
    let mut chunk_index = 0;
    loop {
        // Read a chunk of records.
        let records = read_chunk(&mut reader, args.chunk_size)?;
        if records.is_empty() {
            break;
        }

        // Parallelize counting over the chunk's records.
        let partial_map = records
            .par_iter()
            .fold(
                || FnvHashMap::default(),
                |local_map, rec| merge_kmer_maps(local_map, count_kmers_in_records(std::slice::from_ref(rec), args.k)),
            )
            .reduce(
                || FnvHashMap::default(),
                |global_map, local_map| merge_kmer_maps(global_map, local_map),
            );

        // Build a minimal de Bruijn from this chunk's k-mers, above threshold = 1 (since we might want to do partial merges).
        let partial_dbg = build_debruijn(&partial_map, 1);

        // Serialize partial de Bruijn to disk.
        let chunk_file = args.partial_outdir.join(format!("partial_debruijn_{}.bin", chunk_index));
        let file = File::create(&chunk_file)
            .with_context(|| format!("Failed to create partial de Bruijn file {:?}", chunk_file))?;
        bincode::serialize_into(BufWriter::new(file), &partial_dbg)
            .with_context(|| format!("Failed to write partial de Bruijn to {:?}", chunk_file))?;

        println!("Processed chunk {} with {} records, wrote partial de Bruijn to {:?}", 
            chunk_index, records.len(), chunk_file);
        chunk_index += 1;
    }

    // Merge all partial de Bruijn graphs into one.
    let dir_entries = std::fs::read_dir(&args.partial_outdir)
        .with_context(|| format!("Failed to read partial output directory at {:?}", args.partial_outdir))?;

    let mut final_graph = DeBruijn::new();
    for entry in dir_entries {
        let path = entry?.path();
        if path.file_name().map_or(false, |f| f.to_string_lossy().starts_with("partial_debruijn_")) {
            let partial_graph = read_debruijn_graph(&path)
                .with_context(|| format!("Failed to read partial de Bruijn file {:?}", path))?;
            final_graph.merge(partial_graph);
        }
    }

    // Filter out edges with count below the user threshold (i.e., finalize the graph).
    // We'll do a quick rebuild step for thresholding.
    let mut thresholded_graph = DeBruijn::new();
    for (prefix, edges) in final_graph.adjacency {
        let mut new_edges = FnvHashMap::default();
        for (base, count) in edges {
            if count >= args.threshold {
                *new_edges.entry(base).or_insert(0) += count;
            }
        }
        if !new_edges.is_empty() {
            thresholded_graph.adjacency.insert(prefix, new_edges);
        }
    }

    // Write out the final merged graph.
    write_debruijn_graph(&thresholded_graph, &args.final_output)?;
    println!(
        "Final de Bruijn graph has {} prefix nodes. Written to {:?}.",
        thresholded_graph.adjacency.len(),
        args.final_output
    );

    Ok(())
}
