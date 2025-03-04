## Explanation of the Code and Output:

This program constructs a pangenome graph from multiple FASTA files, each representing a haplotype or genome. The program uses Rayon to parallelize the construction process, making it scalable for large datasets.
Key Components:

    PGNode Struct:
        The PGNode struct represents a node in the pangenome graph. Each node is a k-mer (subsequence of length k).
        The struct implements the Eq, PartialEq, and Hash traits to allow nodes to be used as keys in a HashMap.

    PGGraph Type Alias:
        The PGGraph is a HashMap where each key is a PGNode (representing a k-mer), and the value is a vector of successor PGNodes. This represents directed edges between k-mers.

    build_pangenome_graph Function:
        This function takes a list of FASTA file paths and a k-mer size k. It reads each FASTA file, processes the sequences, and constructs partial k-mer graphs for each file.
        It uses Rayon's parallel iteration (par_iter) to read and process the sequences concurrently, speeding up the construction process.
        For each sequence in a FASTA file, the program:
            Skips sequences shorter than k.
            Extracts consecutive k-mers and their subsequent k-mers to build edges in the graph.
        The function then merges all partial graphs into a single global pangenome graph using parallel reduction (reduce), efficiently combining all the local graphs.

    Main Function:
        The main function defines the paths to the haplotype FASTA files (haplotype1.fasta, haplotype2.fasta, and haplotype3.fasta) and the k-mer size (k = 10).
        It calls build_pangenome_graph to construct the graph and then prints the number of nodes (unique k-mers) in the pangenome graph.
        It also prints a small subset of the graph, showing the k-mers and their successor k-mers (edges).

### Output Breakdown:

    Reading and Processing Files:

Reading file: src/haplotype1.fasta
Reading file: src/haplotype2.fasta
Processing sequence: GATCGATCGATCGATCGTAC
Processing sequence: ATCGATCGATCGATCGATCG
Processing sequence: TACGTAGCTAGCTAGCTAGC
Processing sequence: CGTACGTACGTACGTACGTA
Reading file: src/haplotype3.fasta
Processing sequence: ACGTACGTAGCTAGCTAGCT
Processing sequence: CGTACGATCGTACGTAGCTA

    This part of the output shows the progress as the program reads the input FASTA files and processes each sequence.
    The println! statements inside the map function provide feedback on which file and sequence are currently being processed.

### Constructed Pangenome Graph:

Constructed a pangenome graph with 29 nodes.

    This line indicates that the final pangenome graph contains 29 unique k-mers (nodes). The graph captures the relationships (edges) between these k-mers across the haplotypes provided in the FASTA files.

Small Subset of the Graph:

    Node: CGTAGCTAGC -> ["GTAGCTAGCT", "GTAGCTAGCT"]
    Node: ATCGTACGTA -> ["TCGTACGTAG"]
    Node: GTACGTAGCT -> ["TACGTAGCTA", "TACGTAGCTA"]
    Node: TACGTACGTA -> ["ACGTACGTAC", "ACGTACGTAC"]
    Node: ACGATCGTAC -> ["CGATCGTACG"]

        This section prints a subset of the constructed graph, showing a few k-mers (nodes) and their successor k-mers (edges).
        For example:
            The node CGTAGCTAGC has two outgoing edges to GTAGCTAGCT.
            The node ATCGTACGTA has an edge to TCGTACGTAG.
        This shows how the k-mers are connected based on their overlap in the sequences from the haplotypes.

### Conclusion:

    The program effectively builds a pangenome graph by extracting k-mers from multiple FASTA files (each representing a haplotype or genome).
    The use of Rayon enables parallel processing of the sequences, improving performance when dealing with large datasets.
    The pangenome graph captures the relationships between overlapping k-mers, which is crucial for tasks like genome assembly or comparative genomics.
    The program prints out useful information about the constructed graph, including the total number of unique k-mers and a small subset of the graph with nodes and edges.

This structure forms the basis for more advanced bioinformatics tasks that involve graph-based representations of genomic sequences, allowing for efficient storage and analysis of overlapping sequence data.
