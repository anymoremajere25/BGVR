# Explanation of the Helper Function for Loading a Genomic Sequence
## experiment_2_6_1
The helper function load_genome(path: &str) -> Result<String> is designed to read a genomic sequence from a file and store it in a String for further processing. This function is particularly useful for handling genomic data stored in FASTA format, which consists of a header line (starting with >) followed by nucleotide sequences.
Function Breakdown:

   ### Opening the File:

let file = File::open(path)?;

This line opens the file at the specified path. The ? operator propagates any errors that might occur (e.g., if the file is missing or inaccessible).

## Reading Line by Line Using a Buffered Reader:

let reader = BufReader::new(file);

A BufReader is used to efficiently read the file line by line, reducing memory usage compared to reading the entire file at once.

## *Processing the Genome Data:*

let mut genome = String::new();
for line in reader.lines() {
    let line = line?;
    if !line.starts_with('>') {
        genome.push_str(&line.trim());
    }
}

    The function loops through each line in the file.

    If the line starts with >, it is ignored because it represents a FASTA header, which contains metadata about the sequence rather than nucleotide data.

    Otherwise, the line (trimmed to remove whitespace) is appended to the genome string.

Returning the Loaded Genome Sequence:

    Ok(genome)

    After processing, the entire genome sequence is returned as a single String.

### *FASTA File Source and Explanation*

The genomic sequence file used (complete_seq_100.txt) contains data from NCBI (National Center for Biotechnology Information), specifically complete sequences of SARS-CoV-2 isolates. The FASTA headers in the file indicate that these sequences correspond to different SARS-CoV-2 variants isolated from human samples in the USA in 2025.

*Example of a FASTA header from the file:*

>PV297996.1 Severe acute respiratory syndrome coronavirus 2 isolate SARS-CoV-2/human/USA/CT-DPH-1276957001/2025, complete genome

This indicates:

    PV297996.1: The accession number for this sequence in NCBI.

    Severe acute respiratory syndrome coronavirus 2: The virus (SARS-CoV-2).

    Isolate SARS-CoV-2/human/USA/CT-DPH-1276957001/2025: The specific strain, location, and year of isolation.

    Complete genome: The full-length viral genome is included.

This means the genomic data in complete_seq_100.txt is a full-length SARS-CoV-2 sequence retrieved from NCBI.
Why This Function is Useful

    It simplifies reading large genomic sequences.

    It removes FASTA headers, leaving only nucleotide data for downstream bioinformatics analysis.

    It efficiently handles multi-line genomic sequences, ensuring they are stored as a contiguous string.

  ## Conclusion.

The load_genome function successfully reads and processes a complete genomic sequence from an NCBI FASTA file, stripping out metadata and concatenating the nucleotide sequence into a single string. This processed sequence is then used to construct a suffix array, which is a fundamental data structure for efficient pattern matching and comparative genomics.

### Key Output Insights:

Genome length: 2972496

Suffix array length: 2972496

First 10 entries in suffix array: [2845856, 2311070, 1865468, 2221997, 2073539, 2014076, 2519108, 824517, 2964695, 2281379]

    Genome length (2,972,496 nucleotides): This confirms that the sequence was successfully loaded and processed.

    Suffix array length (2,972,496 indices): Matches the genome length, indicating a complete suffix array was generated.

    First 10 suffix array entries: These indices represent the lexicographically sorted suffixes of the genome, useful for fast substring searches.
