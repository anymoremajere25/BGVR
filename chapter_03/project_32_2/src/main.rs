use mpi::traits::*;
use mpi::environment;
use std::collections::BTreeMap;
use serde::{Serialize, Deserialize};
use bincode;

/// Represents a minimal “partial” suffix array result for a chunk of text.
#[derive(Debug, Clone, Serialize, Deserialize)]
struct PartialSuffixArray {
    rank_id: i32,
    offset: usize,
    sorted_suffixes: Vec<usize>, // Suffix positions relative to offset
}

/// Constructs a naive local suffix array for a chunk of text.
fn build_local_suffix_array(chunk: &str, rank: i32, global_offset: usize) -> PartialSuffixArray {
    let mut suffixes: Vec<usize> = (0..chunk.len()).collect();
    suffixes.sort_by_key(|&pos| &chunk[pos..]);
    PartialSuffixArray {
        rank_id: rank,
        offset: global_offset,
        sorted_suffixes: suffixes,
    }
}

/// Merges partial suffix arrays into a single global array, adjusting offsets.
fn merge_suffix_arrays(chunks: Vec<PartialSuffixArray>, global_text: &str) -> Vec<usize> {
    let mut mapping = BTreeMap::new();
    for c in chunks {
        for &local_pos in &c.sorted_suffixes {
            let global_pos = c.offset + local_pos;
            let suffix_str = &global_text[global_pos..];
            mapping.insert(suffix_str, global_pos);
        }
    }
    mapping.values().copied().collect()
}

fn main() {
    // Initialize the MPI environment
    let _mpi = environment::initialize().expect("Failed to initialize MPI");
    let world = _mpi.world();
    let rank = world.rank();
    let size = world.size();

    // Example text. In real usage, rank 0 might load from a file or read input.
    let text = "GATTACAGATTACACAT";
    let total_len = text.len();

    // Partition text among ranks
    let chunk_size = (total_len + (size as usize) - 1) / (size as usize);
    let start = (rank as usize) * chunk_size;
    let end = ((rank as usize) + 1) * chunk_size;
    let end = end.min(total_len);

    let local_chunk = if start < total_len {
        &text[start..end]
    } else {
        ""
    };

    // Build partial suffix array
    let partial_sa = build_local_suffix_array(local_chunk, rank, start);

    if rank == 0 {
        // Rank 0: collect partial arrays from all ranks into one vector
        let mut partial_arrays = Vec::with_capacity(size as usize);
        partial_arrays.push(partial_sa);

        // Receive partial arrays from ranks 1..(size-1)
        for r in 1..size {
            let (recv_bytes, _status) = world.process_at_rank(r).receive_vec::<u8>();
            let partial: PartialSuffixArray =
                bincode::deserialize(&recv_bytes).expect("Failed to deserialize PartialSuffixArray");
            partial_arrays.push(partial);
        }

        // Merge them into a global suffix array
        let global_sa = merge_suffix_arrays(partial_arrays, text);
        println!("Global Suffix Array: {:?}", global_sa);
    } else {
        // Non-zero ranks serialize PartialSuffixArray and send it to rank 0
        let send_bytes =
            bincode::serialize(&partial_sa).expect("Failed to serialize PartialSuffixArray");
        world.process_at_rank(0).send(&send_bytes[..]);
    }
}