## 4.3 MOTIF DISCOVERY AND REGULATORY ELEMENT IDENTIFICATION
[project_43_2]

This Rust program implements a basic Expectation-Maximization (EM) algorithm for motif discovery in DNA sequences using a Position Weight Matrix (PWM) model. It simulates one EM iteration on an HPC (High-Performance Computing) system.
### 1. What This Code Does

    It initializes a random motif model (PWM) with probabilities for nucleotides (A, C, G, T) at each position.

    It runs an Expectation step (E-step) by scanning DNA sequences to find the best alignment for a motif.

    It performs an approximate Maximization step (M-step) by updating the PWM based on observed nucleotide counts.

    It prints the initial and updated PWM after one iteration.

### 2. Key Components
(a) MotifModel Struct

    Stores the Position Weight Matrix (PWM) as a vector of arrays (pwm).

    Implements:

        new_random(motif_length): Initializes a random PWM.

        score_kmer(kmer): Calculates the probability of a k-mer using the PWM.

        update_from_counts(partial_counts): Updates PWM based on new observations.

(b) EM Iteration

    em_iteration(model, sequences): Simulates the Expectation step (E-step) by:

        Scanning each sequence for the best motif alignment.

        Counting nucleotide occurrences in the best-aligned motifs.

        Returning partial counts, which will later update the PWM.

(c) Main Function

    Initializes a random 4-length motif model.

    Processes three DNA sequences.

    Performs one EM iteration, updating the PWM.

    Prints the initial and updated PWM.

### 3. Output Explanation
Initial PWM (Randomly Initialized)

Initial motif model: 
[
 [0.056, 0.139, 0.350, 0.454], 
 [0.011, 0.394, 0.226, 0.367], 
 [0.177, 0.475, 0.299, 0.047], 
 [0.308, 0.174, 0.244, 0.273]
]

    Each row represents a position in the motif.

    Each column represents A, C, G, T probabilities (e.g., in the first position, A=5.6%, C=13.9%, G=35%, T=45.4%).

Updated PWM (After One EM Iteration)

Updated motif model: 
[
 [0.0, 0.0, 0.333, 0.666], 
 [0.0, 0.0, 0.333, 0.666], 
 [0.666, 0.0, 0.333, 0.0], 
 [0.0, 0.666, 0.333, 0.0]
]

    After observing the sequences, the PWM is updated:

        The first two positions favor G (33%) and T (66%).

        The third position favors A (66%) and G (33%).

        The last position favors C (66%) and G (33%).

### 4. Summary

    The program finds motifs in DNA sequences using probability-based pattern recognition.

    It applies one step of the EM algorithm to refine the motif model.

    Future iterations would further refine PWM for better motif detection.
