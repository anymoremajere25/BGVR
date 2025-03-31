## 4.6 SINGLE-CELL FUNCTIONAL GENOMICS
### Explanation of the CSR Matrix Multiplication in Rust [experiment_4_6_2]

This Rust program implements parallel multiplication of a Compressed Sparse Row (CSR) matrix with a vector using Rayon for parallelism. The output of the program is the result of this multiplication.
### 1. What is a CSR Matrix?

A Compressed Sparse Row (CSR) matrix is an efficient way to store sparse matrices (matrices with mostly zero values) to save memory and improve performance. It consists of:

    values: Stores non-zero values in row-major order.
    col_indices: Stores the column indices corresponding to the non-zero values.
    row_ptrs: Marks where each row starts and ends in the values and col_indices arrays.


#### 2. Step-by-Step Breakdown
Step 1: Define the CSR Matrix

let csr = CsrMatrix {
    values: vec![5.0, 2.0, 3.0],  // Non-zero elements
    col_indices: vec![0, 2, 3],   // Column indices for each non-zero value
    row_ptrs: vec![0, 1, 2, 3],   // Pointers to the start of each row
    nrows: 3,                     // Number of rows
    ncols: 4,                     // Number of columns
};

This matrix represents:
[5.00.00.00.00.00.02.00.00.00.00.03.0]
​5.00.00.0​0.00.00.0​0.02.00.0​0.00.03.0​
​
Step 2: Define the Input Vector

let vec = vec![2.0, 0.0, 1.0, 1.5];

This is the vector that we multiply with the CSR matrix.
Step 3: Perform Parallel Multiplication

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
        let mut guard = result.lock().unwrap();
        guard[row] = sum;
    });

    Arc::try_unwrap(result).unwrap().into_inner().unwrap()
}

How It Works

    Arc<Mutex<Vec<f64>>>
        We create a shared vector (result) to store the output.
        Arc enables shared memory across threads.
        Mutex ensures safe concurrent updates.

    Parallel Row Processing
        into_par_iter() enables parallel execution.
        Each row computes its result independently.

    Iterate Over Non-Zero Values
        start_ptr and end_ptr find where the row's non-zero values are stored.
        Multiply each non-zero matrix element with the corresponding vector element.
        Accumulate the sum for the row.

    Update Output Vector Safely
        The final result for the row is written inside a locked block to prevent race conditions.

### Step 4: Print the Result

let product = csr.mul_vector_parallel(&vec);
println!("CSR * vector = {:?}", product);

The output is:

CSR * vector = [10.0, 2.0, 4.5]

Calculation Breakdown

For each row:

    Row 0: 5.0×2.0=10.05.0×2.0=10.0
    Row 1: 2.0×1.0=2.02.0×1.0=2.0
    Row 2: 3.0×1.5=4.53.0×1.5=4.5

### 5. Summary

✔ CSR format reduces memory usage by storing only non-zero values.
✔ Parallelism speeds up matrix-vector multiplication using Rayon.
✔ Mutex ensures safe concurrent updates to avoid race conditions.
⚠ Locking overhead can be optimized → Lock-free approaches (e.g., per-thread accumulation) may improve performance.

