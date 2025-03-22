use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;

/// A minimal HyperLogLog implementation in Rust.
#[derive(Clone)]
pub struct HyperLogLog {
    p: u32,             // precision, such that m = 2^p
    registers: Vec<u8>, // each register holds the max run-of-zeros for a bucket
    alpha: f64,         // bias-correcting constant
}

impl HyperLogLog {
    /// Create a new HyperLogLog structure with precision p.
    /// This yields m = 2^p registers.
    pub fn new(p: u32) -> Self {
        let m = 1 << p;
        // The alpha constant depends on m. 
        let alpha = match m {
            16 => 0.673,
            32 => 0.697,
            64 => 0.709,
            _  => 0.7213 / (1.0 + 1.079 / m as f64),
        };
        HyperLogLog {
            p,
            registers: vec![0; m],
            alpha,
        }
    }

    /// Convenience constructor that builds a HyperLogLog from an iterator of items.
    /// Each item in the iterator is added to the HyperLogLog.
    pub fn from_iter<I, T>(p: u32, iter: I) -> Self
    where
        I: IntoIterator<Item = T>,
        T: Hash,
    {
        let mut hll = Self::new(p);
        for item in iter {
            hll.add(&item);
        }
        hll
    }

    /// Hash an item to produce a 64-bit value.
    fn hash<T: Hash>(item: &T) -> u64 {
        let mut hasher = DefaultHasher::new();
        item.hash(&mut hasher);
        hasher.finish()
    }

    /// Add an element to the HyperLogLog, updating the relevant register.
    pub fn add<T: Hash>(&mut self, item: &T) {
        let hash_val = Self::hash(item);
        // Extract top p bits for the bucket.
        let bucket = (hash_val >> (64 - self.p)) as usize;
        // Discard top p bits and count leading zeros in the remainder.
        let w = hash_val << self.p;
        let leading_zeros = (w.leading_zeros() + 1) as u8;
        // Store maximum run of zeros observed for this bucket.
        self.registers[bucket] = self.registers[bucket].max(leading_zeros);
    }

    /// Merge another HyperLogLog into this one (they must have the same precision p).
    /// This allows you to combine counts from multiple data streams.
    pub fn merge(&mut self, other: &Self) {
        assert_eq!(
            self.p, other.p,
            "Cannot merge HyperLogLogs with different precisions."
        );
        for (a, b) in self.registers.iter_mut().zip(other.registers.iter()) {
            *a = (*a).max(*b);
        }
    }

    /// Estimate the cardinality of the set using the HyperLogLog formula.
    pub fn estimate(&self) -> f64 {
        let m = 1 << self.p;
        let sum: f64 = self.registers
            .iter()
            .map(|&val| 2_f64.powi(-(val as i32)))
            .sum();
        // alpha * m^2 / sum
        let raw_est = self.alpha * (m as f64) * (m as f64) / sum;
        raw_est
    }
}

/// A small demonstration of HyperLogLog usage with sample data.
fn main() {
    // Example 1: Simple list of integer values
    let int_values: Vec<u32> = (0..10_000).collect();
    let mut hll_int = HyperLogLog::from_iter(10, int_values.clone());
    let estimate_int = hll_int.estimate();
    println!("Actual integer count: {}", int_values.len());
    println!("Estimated integer count: {:.2}", estimate_int);

    // Example 2: Simple list of string values (with duplicates)
    let strings = vec!["apple", "banana", "cherry", "banana", "date", "apple", "elderberry"];
    let mut hll_str = HyperLogLog::from_iter(4, strings.iter());
    let estimate_str = hll_str.estimate();
    // Calculate the actual unique values
    let actual_unique_str = strings.iter().collect::<std::collections::HashSet<_>>().len();
    println!("\nActual string unique count: {}", actual_unique_str);
    println!("Estimated string count: {:.2}", estimate_str);

    // Example 3: Demonstrate merging
    // We'll split int_values into two parts, create two separate HyperLogLogs, then merge.
    let (part1, part2) = int_values.split_at(5000);
    let hll_part1 = HyperLogLog::from_iter(10, part1.iter());
    let hll_part2 = HyperLogLog::from_iter(10, part2.iter());

    // Merge them into a new structure
    let mut hll_merged = hll_part1.clone();
    hll_merged.merge(&hll_part2);

    let merged_estimate = hll_merged.estimate();
    println!("\nMerging two HyperLogLogs each containing half of the integer range:");
    println!("Merged estimate of unique integers: {:.2}", merged_estimate);
    println!("Actual unique count (0..10000): {}", 10_000);
}
