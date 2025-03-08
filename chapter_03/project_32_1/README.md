## 3.2. Sequence Data Structures and Strings Algorithms##
### project__32_1

This code demonstrates a parallelized approach for detecting a specific genomic pattern (such as a short sequence) within large genomic datasets using MPI and the Knuth-Morris-Pratt (KMP) algorithm. The genomic reference or read data is divided among multiple processes, ensuring each node handles a portion of the text along with necessary overlaps. This prevents missing matches that span across chunk boundaries. The approach scales well across numerous processes, making it ideal for high-performance computing (HPC) environments where genomic references can reach gigabytes or terabytes in size.

In this implementation, Rank 0 loads the genomic data from a file, calculates chunk boundaries, and distributes the appropriate sections (including any required overlaps) to other ranks. Each rank performs a local KMP search on its assigned chunk and reports partial match indices back to Rank 0. Since these match positions are relative to individual chunks, Rank 0 adjusts them accordingly when assembling the final list of matches. The partial results are then merged and sorted to generate a global list of pattern occurrences.
Project Structure

project_32_1/

    Cargo.toml – Configuration file specifying dependencies
   project_32_1/src/
        main.rs – Rust implementation of the algorithm
        big_text_example.txt – Sample dataset (text file)
        generate_big_text_example.txt.ipynb – Python script to generate the dataset
        output.txt – Output file storing the results

## How to Run

To execute the program in PowerShell, run:

cargo run > output.txt

(This runs main.rs and saves the output to output.txt.)
Dependencies

rand = "0.9.0"  
rayon = "1.10.0"

### Output Explanation

The program’s output lists the global positions in the text where the specified pattern "ABABABC" appears. These positions are determined using a parallelized string-matching approach, leveraging the KMP algorithm and MPI for distributed execution.
1. Global Matches for Pattern "ABABABC"

    The pattern "ABABABC" is detected multiple times throughout the text.
    The output lists the starting indices where the pattern appears.
    There are 200 matches, with each occurrence approximately 5000 characters apart, suggesting a repeating structure within the text.

2. Text Processing and Parallelization

    The program utilizes MPI to divide the text into chunks, distributing the workload across multiple processes.
    Rank 0 coordinates the data distribution, sending chunks to different ranks for independent processing.
    Each rank performs a KMP search on its assigned chunk and returns local match positions, which are later merged by Rank 0.

3. Pattern Matching

    The KMP algorithm efficiently locates occurrences of "ABABABC" by precomputing a prefix table, which allows it to skip unnecessary comparisons.
    When a match is found, its starting index is recorded.
    These local matches are consolidated into a final sorted list of global pattern occurrences.

### Conclusion

The results indicate that the pattern "ABABABC" appears consistently within the dataset at positions such as 4993, 9993, 14993, and so on, with intervals of approximately 5000 characters.
Key Takeaways:

    The program efficiently searches large texts by leveraging MPI for parallel computing.
    The KMP algorithm ensures optimized pattern matching, making it suitable for large-scale datasets.
    The final output represents the global match positions after merging results from all processes.
    This approach is highly scalable and can be applied to text mining and other domains requiring fast pattern detection in massive datasets.
