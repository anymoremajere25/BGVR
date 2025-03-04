## Explanation of the Output

### De Bruijn Graph Construction:

Constructed De Bruijn graph with 19 nodes.

    The program builds a De Bruijn graph using a k-mer size of 21 (as specified in the code).
    The graph is constructed by processing all the sequences in parallel using the Rayon library, which allows multi-threading to speed up the computation.
    Each sequence is processed to find all possible k-mers of length 21. For each k-mer, the program records the subsequent k-mer (i.e., the k-mer that starts one base after the current k-mer). This creates an edge from one k-mer to another.
    After processing the sequences, the program merges local HashMaps from each thread to produce a global De Bruijn graph. The graph consists of 19 nodes (unique k-mers), which means there are 19 distinct k-mers found in the sequences.

### Matrix Dimensions (nalgebra & ndarray):

    nalgebra matrix: 5 x 5
    ndarray shape: 5 x 5

        These two lines display the dimensions of the matrices created using the nalgebra and ndarray crates.
            The nalgebra matrix is a 5x5 matrix filled with the value 1.0. This is created using DMatrix::<f32>::from_element(5, 5, 1.0).
            The ndarray matrix is also 5x5, created using Array2::<f32>::ones((5, 5)).
            These matrices are just examples to demonstrate how matrices can be created using the nalgebra and ndarray libraries. They do not directly relate to the construction of the De Bruijn graph but are included for illustration purposes (likely for later use in high-performance computations).

### Conclusion:

    De Bruijn Graph Construction:
        The program constructs a De Bruijn graph from the sequences in the FASTA file. A k-mer size of 21 is chosen, and for each sequence, k-mers and their subsequent k-mers (edges) are identified. The graph contains 19 unique nodes (k-mers), meaning that there are 19 distinct sequences of length 21 found in the input sequences.
        The use of Rayon for parallel processing allows the program to scale efficiently by splitting the work of constructing the graph across multiple threads, which is particularly useful when dealing with a large number of sequences.

  ### Matrix Operations:
        The program also demonstrates the creation of 5x5 matrices using the nalgebra and ndarray libraries. These matrices are initialized with 1.0 and 1 values, respectively. This part of the program is more of an example for future high-performance computing tasks rather than directly related to the graph construction.

In summary, the program efficiently constructs a De Bruijn graph from DNA sequences in parallel, demonstrating scalability using Rayon for multi-threading. It also introduces basic matrix operations, which could be useful in further numerical or computational tasks.
