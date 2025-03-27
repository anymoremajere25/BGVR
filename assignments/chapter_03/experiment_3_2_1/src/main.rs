use mpi::environment::Universe;
use mpi::traits::*;
use std::error::Error;
use std::fs;

/// Builds the prefix table (failure function) for the KMP algorithm.
fn build_prefix_table(pattern: &str) -> Vec<usize> {
    let pchars: Vec<char> = pattern.chars().collect();
    let mut prefix = vec![0; pchars.len()];
    let mut j = 0; // length of the previous longest prefix suffix

    // Loop calculates prefix[i] for i in [1..pattern.len())
    for i in 1..pchars.len() {
        while j > 0 && pchars[i] != pchars[j] {
            j = prefix[j - 1];
        }
        if pchars[i] == pchars[j] {
            j += 1;
            prefix[i] = j;
        }
    }
    prefix
}

/// KMP search over a local text chunk, returning match positions (relative to chunk start).
fn kmp_search_local(pattern: &str, text_chunk: &str) -> Vec<usize> {
    let prefix = build_prefix_table(pattern);
    let pchars: Vec<char> = pattern.chars().collect();
    let tchars: Vec<char> = text_chunk.chars().collect();
    let mut matches = Vec::new();

    let mut j = 0; // index into pattern
    for i in 0..tchars.len() {
        while j > 0 && tchars[i] != pchars[j] {
            j = prefix[j - 1];
        }
        if tchars[i] == pchars[j] {
            j += 1;
        }
        if j == pchars.len() {
            // Found a match ending at i => match starts at i+1 - pchars.len()
            matches.push(i + 1 - pchars.len());
            j = prefix[j - 1];
        }
    }
    matches
}

/// Compute chunk boundaries for each rank using simple partitioning.
fn compute_chunk_boundaries(text_len: usize, size: usize) -> Vec<(usize, usize)> {
    // chunk_size = ceil(text_len / size)
    let chunk_size = (text_len + size - 1) / size;
    let mut boundaries = Vec::new();
    for rank in 0..size {
        let start = rank * chunk_size;
        let mut end = start + chunk_size;
        if end > text_len {
            end = text_len;
        }
        boundaries.push((start, end));
    }
    boundaries
}

fn main() -> Result<(), Box<dyn Error>> {
    let universe: Universe = mpi::initialize().expect("Failed to initialize MPI");
    let world = universe.world();
    let rank = world.rank();
    let size = world.size() as usize;

    // Example pattern
    let pattern = "ABABABC";

    // Rank 0 loads text from file or fallback
    let text_path = "big_text_example.txt";
    let text = if rank == 0 {
        match fs::read_to_string(text_path) {
            Ok(s) => s,
            Err(_) => {
                // fallback if file not found
                "ABABABCAABABABABCABABABCA".to_string()
            }
        }
    } else {
        String::new()
    };

    // Rank 0 calculates chunk boundaries
    let text_len = if rank == 0 { text.len() } else { 0 };

    // Broadcast text length to all ranks
    let mut buf_len = vec![text_len as u64];
    world.process_at_rank(0).broadcast_into(&mut buf_len[..]);
    let global_text_len = buf_len[0] as usize;

    // Rank 0 => compute chunk boundaries
    let chunk_boundaries = if rank == 0 {
        compute_chunk_boundaries(global_text_len, size)
    } else {
        Vec::new()
    };

    // Distribute chunks:
    //  1) Rank 0: for each rank r, sends [start, end] as a Vec<u64>, then chunk bytes
    //  2) Ranks > 0: receive the boundary vector, then the chunk bytes
    let (my_start, _my_end, local_text) = if rank == 0 {
        let mut local = (0, 0, String::new());
        for r in 0..size {
            let (st, en) = chunk_boundaries[r];
            let chunk_str = &text[st..en];
            let chunk_bytes = chunk_str.as_bytes().to_vec();

            if r as i32 == 0 {
                // rank 0 uses it locally
                local = (st, en, chunk_str.to_string());
            } else {
                // send start, end as a dynamic Vec<u64>
                let offset_vec = vec![st as u64, en as u64];
                world.process_at_rank(r as i32).send(&offset_vec[..]);
                // send chunk bytes
                world.process_at_rank(r as i32).send(&chunk_bytes[..]);
            }
        }
        local
    } else {
        // Non-root ranks: receive the offset vec, then chunk
        let (offset_vec, _status) = world.process_at_rank(0).receive_vec::<u64>();
        let st = offset_vec[0] as usize;
        let en = offset_vec[1] as usize;
        let (chunk_bytes, _status) = world.process_at_rank(0).receive_vec::<u8>();
        (st, en, String::from_utf8_lossy(&chunk_bytes).to_string())
    };

    // Perform KMP on local_text
    let local_matches = kmp_search_local(pattern, &local_text);

    // Each rank sends partial matches back to rank 0:
    // We'll do a 2-step:
    //   1) send start offset in a Vec<u64> of length 1
    //   2) send local matches in a Vec<u64>
    if rank != 0 {
        // 1) offset
        let offset_vec = vec![my_start as u64];
        world.process_at_rank(0).send(&offset_vec[..]);
        // 2) local matches
        let matches_vec = local_matches.iter().map(|&x| x as u64).collect::<Vec<_>>();
        world.process_at_rank(0).send(&matches_vec[..]);
    } else {
        // rank 0 merges
        let mut global_matches = Vec::new();
        // Add local rank 0's matches
        for &pos in &local_matches {
            global_matches.push(pos + my_start);
        }
        // Receive from other ranks
        for r in 1..size {
            // step 1) offset
            let (offset_vec, _) = world.process_at_rank(r as i32).receive_vec::<u64>();
            let st = offset_vec[0] as usize;
            // step 2) partial matches
            let (matches_vec, _) = world.process_at_rank(r as i32).receive_vec::<u64>();
            for pm in matches_vec {
                global_matches.push(pm as usize + st);
            }
        }

        global_matches.sort_unstable();
        println!("Global matches for pattern '{p}' => {m:?}",
                 p = pattern, m = global_matches);
    }

    Ok(())
}