## EXPLANATION /OVERVIEW OF project_31

This Rust program demonstrates the use of a Bloom Filter to efficiently test membership in a set of DNA sequences. It uses parallelism to speed up the insertion of sequences into the filter, and it handles random sequence generation and testing.
### 1. Bloom Filter Implementation:
### 1.1 Bloom Filter Struct:

pub struct BloomFilter {
    bits: Vec<AtomicBool>,  // Vector of AtomicBool to store the bits in the filter
    k: usize,               // Number of hash functions (or seeds)
    m: usize,               // Size of the bit array
}

    bits: A vector of AtomicBool is used to store bits of the Bloom Filter. AtomicBool ensures thread-safe updates to the bits without using locks.
    k: The number of hash functions used to hash items into the Bloom Filter. Each item will be hashed k times with different seeds.
    m: The size of the bit array. The Bloom Filter has m bits, and the range of indices for the bits is 0..m.

### 1.2 Creating a New Bloom Filter:

pub fn new(m: usize, k: usize) -> Self {
    let bits = (0..m).map(|_| AtomicBool::new(false)).collect();
    Self { bits, k, m }
}

    This function creates a new Bloom Filter with m bits, where each bit is initially set to false.
    The filter also receives k hash functions (or seed values).

### 1.3 Inserting an Item:

pub fn insert<T: Hash + Sync>(&self, item: &T) {
    (0..self.k).into_par_iter().for_each(|seed| {
        let pos = self.hash_with_seed(item, seed) % self.m;
        self.bits[pos].store(true, Ordering::Relaxed);
    });
}

    Parallel Insertion: This function hashes the item k times using different seeds. Each hash generates a bit position. It then sets the bit at that position to true.
    Parallelism: The function uses into_par_iter() from the rayon crate to parallelize the hashing and setting of bits across k iterations. This speeds up the process when dealing with a large number of items.

### 1.4 Checking for Membership:

pub fn contains<T: Hash + Sync>(&self, item: &T) -> bool {
    (0..self.k).into_par_iter().all(|seed| {
        let pos = self.hash_with_seed(item, seed) % self.m;
        self.bits[pos].load(Ordering::Relaxed)
    })
}

    This function checks if an item might be in the Bloom Filter. It hashes the item k times and checks if the bits at the resulting positions are all true. If any bit is false, it returns false because the item is definitely not in the filter.
    Again, the function uses parallelism (into_par_iter()) to check the positions for each hash.

### 1.5 Hashing Function with Seed:

### [inline]
fn hash_with_seed<T: Hash>(&self, item: &T, seed: usize) -> usize {
    let mut hasher = DefaultHasher::new();
    item.hash(&mut hasher);   // Hash the item
    seed.hash(&mut hasher);   // Mix in the seed
    hasher.finish() as usize
}

    This function hashes the item combined with a seed. The resulting hash is used to determine the position in the Bloom Filter's bit array.

### 2. Random DNA Sequence Generation:
2.1 Generating a Random DNA Sequence:

fn generate_random_dna(len: usize) -> String {
    let mut r = rng();
    let bases = ['A', 'C', 'G', 'T'];
    (0..len).map(|_| *bases.choose(&mut r).unwrap()).collect()
}

    This function generates a random DNA sequence of length len. It randomly picks a base ('A', 'C', 'G', or 'T') for each position in the sequence. It uses rng() from the rand crate to get a random number generator and selects bases from the bases array.

### 2.2 Generating Multiple DNA Sequences:

fn generate_synthetic_genomic_data(num_sequences: usize, seq_length: usize) -> Vec<String> {
    (0..num_sequences)
        .map(|_| generate_random_dna(seq_length))
        .collect()
}

    This function generates multiple DNA sequences (num_sequences), each of length seq_length.

### 3. Main Function:
### 3.1 Creating the Bloom Filter:

let bloom_size = 50_000;
let num_hashes = 5;
let bloom = BloomFilter::new(bloom_size, num_hashes);

    A Bloom Filter is created with a bit array size of 50,000 bits and 5 hash functions.

### 3.2 Generating DNA Data:

let num_sequences = 5_000;
let seq_length = 20;
let dataset = generate_synthetic_genomic_data(num_sequences, seq_length);

    5,000 DNA sequences of length 20 are generated using the generate_synthetic_genomic_data function.

3.3 Inserting Sequences into the Bloom Filter:

dataset.par_iter().for_each(|seq| {
    bloom.insert(seq);
});

    The DNA sequences are inserted into the Bloom Filter in parallel using par_iter() from the rayon crate.

### 3.4 Testing Membership:

let test_sequences = vec![
    dataset[0].clone(),
    dataset[1].clone(),
    generate_random_dna(seq_length),
    generate_random_dna(seq_length),
];

for seq in &test_sequences {
    let result = bloom.contains(seq);
    println!("Sequence: {} => in Bloom Filter? {}", seq, result);
}

    The program tests membership for 4 sequences:
        Two sequences that are known to be in the dataset (dataset[0] and dataset[1]).
        Two randomly generated sequences that are likely not in the dataset.
    The contains() function checks if each sequence might be in the Bloom Filter, and the results are printed.

### 4. Output Explanation:

Sequence: TTAAGTAGCGGCTGTTATGT => in Bloom Filter? true
Sequence: CCATCCGCCGTCTTACATAC => in Bloom Filter? true
Sequence: GCAATAATCACCATCCACTA => in Bloom Filter? false
Sequence: CTACGGCCTTAGCAATGAGC => in Bloom Filter? false

   ### First two sequences (TTAAGTAGCGGCTGTTATGT and CCATCCGCCGTCTTACATAC):
        These are sequences that were inserted into the Bloom Filter earlier. They might be in the filter, and since the Bloom Filter is designed to minimize false negatives, they are returned as true (though false positives can still occur in a Bloom Filter, it's unlikely here).

   ### Last two sequences (GCAATAATCACCATCCACTA and CTACGGCCTTAGCAATGAGC):
        These are random sequences, and based on the nature of Bloom Filters, these sequences are definitely not in the filter, so we expect the result to be false. However, since Bloom Filters can produce false positives, we might occasionally get a true, but in this case, both return false.

### Conclusion:

The Bloom Filter is successfully used to test whether a sequence is part of a large set of DNA sequences. The output demonstrates how Bloom Filters can give fast membership tests with a small chance of false positives, but no false negatives. The program also uses parallelism to efficiently handle a large number of sequences.

