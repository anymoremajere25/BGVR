use rayon::prelude::*;
use std::sync::{Arc, Mutex};

/// Basic representation of a CSR (Compressed Sparse Row) matrix.
pub struct CsrMatrix {
    // The non-zero values in row-major order.
    pub values: Vec<f64>,
    // The corresponding column indices for each non-zero value.
    pub col_indices: Vec<usize>,
    // row_ptrs[i] points to the start of row i in values/col_indices.
    pub row_ptrs: Vec<usize>,
    pub nrows: usize,
    pub ncols: usize,
}

impl CsrMatrix {
    /// Perform parallel CSR * vector multiplication using Rust concurrency.
    /// We store the partial sums in a shared vector protected by a Mutex.
    /// For HPC-scale usage, a lock-free approach or local accumulation can reduce contention.
    pub fn mul_vector_parallel(&self, vec: &[f64]) -> Vec<f64> {
        assert_eq!(vec.len(), self.ncols);
        let result = Arc::new(Mutex::new(vec![0.0; self.nrows]));

        (0..self.nrows).into_par_iter().for_each(|row| {
            let start_ptr = self.row_ptrs[row];
            let end_ptr = self.row_ptrs[row + 1];
            let mut sum = 0.0;
            for idx in start_ptr..end_ptr {
                let col_idx = self.col_indices[idx];
                let val = self.values[idx];
                sum += val * vec[col_idx];
            }
            // Acquire lock for final update
            let mut guard = result.lock().unwrap();
            guard[row] = sum;
        });

        Arc::try_unwrap(result).unwrap().into_inner().unwrap()
    }
}

fn main() {
    // Example: 3x4 CSR matrix with minimal non-zero values.
    let csr = CsrMatrix {
        values: vec![5.0, 2.0, 3.0],
        col_indices: vec![0, 2, 3],
        row_ptrs: vec![0, 1, 2, 3],
        nrows: 3,
        ncols: 4,
    };

    // A dense vector for multiplication.
    let vec = vec![2.0, 0.0, 1.0, 1.5];

    // Perform the parallel multiplication.
    let product = csr.mul_vector_parallel(&vec);
    println!("CSR * vector = {:?}", product);
}
