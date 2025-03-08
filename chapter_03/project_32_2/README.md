# 3.2. Sequence Data Structures and Strings Algorithms

## Explanation of MPI-Based Suffix Array Code [project_32_2]

This Rust program parallelizes the construction of a suffix array using MPI (Message Passing Interface). Here's a breakdown of the code:
### 1. What is a Suffix Array?

A suffix array is an array of integers representing the starting positions of all suffixes of a given string, sorted in lexicographical order.

#### For example, given the string "banana", the suffixes are:
Index	Suffix
5	"a"
3	"ana"
1	"anana"
0	"banana"
4	"na"
2	"nana"

 The suffix array would be: [5, 3, 1, 0, 4, 2].
#### 2. Overview of the Program

    The program uses MPI to divide a string across multiple processors.
    Each processor computes a local suffix array for its chunk of the text.
    The partial suffix arrays are gathered and merged into a global suffix array.

#### 3. Key Components of the Code
3.1. Dependencies in Cargo.toml

[dependencies]
mpi = "0.8.0"
serde = { version = "1.0", features = ["derive"] }
bincode = "1.3"

    mpi: Enables parallel computing using the MPI library.
    serde: Serializes and deserializes Rust structs.
    bincode: A binary serialization format used for sending data between processes.

3.2. PartialSuffixArray Struct

#[derive(Debug, Clone, Serialize, Deserialize)]
struct PartialSuffixArray {
    rank_id: i32,
    offset: usize,
    sorted_suffixes: Vec<usize>,
}

    Stores the local suffix array computed by each process.
    rank_id: The MPI rank (ID) of the process computing this array.
    offset: The starting index of the text chunk in the global string.
    sorted_suffixes: A vector of suffix indices relative to the chunk.

3.3. build_local_suffix_array Function

fn build_local_suffix_array(chunk: &str, rank: i32, global_offset: usize) -> PartialSuffixArray {
    let mut suffixes: Vec<usize> = (0..chunk.len()).collect();
    suffixes.sort_by_key(|&pos| &chunk[pos..]);
    PartialSuffixArray {
        rank_id: rank,
        offset: global_offset,
        sorted_suffixes: suffixes,
    }
}

    Takes a substring chunk and sorts suffix positions relative to the chunk.
    Uses sort_by_key to sort suffixes based on their lexicographical order.
    Returns a PartialSuffixArray with the sorted indices.

3.4. merge_suffix_arrays Function

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

    Merges all partial suffix arrays into a global suffix array.
    Uses a BTreeMap to sort suffixes in lexicographical order globally.
    Returns a vector of global suffix positions.

3.5. MPI Execution in main()

fn main() {
    let _mpi = environment::initialize().expect("Failed to initialize MPI");
    let world = _mpi.world();
    let rank = world.rank();
    let size = world.size();

    Initializes MPI.
    Retrieves the rank (ID) and total number of processes.

3.5.1. Partitioning the Text

let text = "GATTACAGATTACACAT";
let total_len = text.len();

let chunk_size = (total_len + (size as usize) - 1) / (size as usize);
let start = (rank as usize) * chunk_size;
let end = ((rank as usize) + 1) * chunk_size;
let end = end.min(total_len);

let local_chunk = if start < total_len {
    &text[start..end]
} else {
    ""
};

    The text is divided into equal-sized chunks for each process.
    Each process extracts its own chunk.

3.5.2. Compute Local Suffix Array

let partial_sa = build_local_suffix_array(local_chunk, rank, start);

    Each process computes a local suffix array.

3.5.3. Collecting and Merging (Rank 0)

if rank == 0 {
    let mut partial_arrays = Vec::with_capacity(size as usize);
    partial_arrays.push(partial_sa);

    for r in 1..size {
        let (recv_bytes, _status) = world.process_at_rank(r).receive_vec::<u8>();
        let partial: PartialSuffixArray = bincode::deserialize(&recv_bytes).expect("Failed to deserialize PartialSuffixArray");
        partial_arrays.push(partial);
    }

    let global_sa = merge_suffix_arrays(partial_arrays, text);
    println!("Global Suffix Array: {:?}", global_sa);
}

    Rank 0 collects suffix arrays from all processes.
    Deserializes them using bincode.
    Merges them into a global suffix array.

3.5.4. Sending Data to Rank 0 (Other Ranks)

else {
    let send_bytes = bincode::serialize(&partial_sa).expect("Failed to serialize PartialSuffixArray");
    world.process_at_rank(0).send(&send_bytes[..]);
}

    Other ranks serialize their local suffix array and send it to rank 0.

#### 4. Output Explanation

Global Suffix Array: [11, 4, 13, 6, 15, 8, 1, 12, 5, 14, 7, 0, 16, 10, 3, 9, 2]

    The output represents sorted suffix positions in the original text.
    Each number is an index in "GATTACAGATTACACAT", sorted in lexicographical order.

#### 5. Key Takeaways

✅ Parallelized with MPI: Each process computes a suffix array for its chunk.
✅ Efficient Communication: Uses bincode for serialization and MPI for sending/receiving data.
✅ Merging with BTreeMap: Ensures correct lexicographical order of suffixes.
✅ Works for Large Inputs: Can handle larger texts with multiple processes.
#### 6. Next Steps

    Test with larger inputs: Modify text with a longer sequence.
    Increase the number of processes:

mpirun -np 4 cargo run --release

Benchmark performance:

time mpirun -np 4 cargo run --release

Optimize sorting: Implement a suffix tree or radix sort for better efficiency.
