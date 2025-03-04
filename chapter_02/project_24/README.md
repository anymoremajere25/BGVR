## Introduction to Rust Programming Language (project_24)

In bioinformatics workflows, it's common to handle files representing various stages of an analysis pipeline, such as FASTQ files containing raw sequencing reads, BAM files showing read alignments to a reference genome, and VCF files listing detected variants. To demonstrate a unified approach, the following Rust code processes each file format in parallel to search for a specific DNA motif (e.g., “GATTACA”). While real-world bioinformatics projects often require more complex operations, such as analyzing CIGAR alignments in BAM files or genotype fields in VCF files, this example highlights how Rust can efficiently handle these formats, leveraging parallel processing with Rayon.

The code defines a helper function, count_occurrences, that identifies overlapping occurrences of a target motif within a string. Three specialized functions—process_fastq, process_bam, and process_vcf—are provided, each designed to process a different file format using community libraries like bio (for FASTQ files) and rust-htslib (for SAM/BAM and VCF files). Within each function, Rust reads the input records, extracts relevant sequences or alleles into a Vec, and utilizes par_iter() from Rayon to distribute the motif-counting workload across available CPU cores. The main function coordinates these tasks by specifying the motif (e.g., “GATTACA”), applying the respective function to each file type, and outputting the total motif matches found. This approach separates the reading logic from the parallel processing, maintains memory safety, and allows easy scaling for larger datasets or more advanced analyses, ensuring both performance and maintainability.
## Files Included:

    main.rs (Rust script)
    main.nf (Nextflow script)
    example.fastq.rara (compressed FASTQ file)
    example.bam (BAM file)
    example.vcf (VCF file)
    Cargo.toml (Cargo configuration file)
    output.txt (Output file)

## How to Run:

Execute the following command in a WSL terminal to run the analysis and save the output:

nextflow run main.nf | tee output.txt

Dependencies:

    bio = "2.0.3"
    rayon = "1.10.0"
    rust-htslib = "0.49.0"

Explanation of the Output:

    FASTQ: Found motif 'GATTACA' 0 times.
        The program processed the example.fastq file, searching for occurrences of "GATTACA". Since the count is 0, it means the motif was not found in any of the sequences. This could be due to the motif being absent in the dataset, sequences being too short, or mutations preventing an exact match.

    BAM: Found motif 'GATTACA' 230 times.
        The program analyzed example.bam, a binary alignment file containing mapped sequencing reads. "GATTACA" was detected 230 times across the aligned reads, indicating that the motif appears frequently in the mapped sequencing data. This could suggest biological significance, such as a conserved sequence or a potential sequencing artifact.

    VCF: Found motif 'GATTACA' 0 times.
        The program processed example.vcf, which contains variant data. Since the motif was not found in any reference or alternative alleles, it suggests that the dataset does not include variations that introduce this motif, or the motif does not overlap with variant positions in the VCF file.

## Conclusion:

The program effectively analyzed sequencing and variant data from different file formats (FASTQ, BAM, and VCF) to identify occurrences of the "GATTACA" motif. Although the motif was absent in the raw sequencing reads (FASTQ) and variant calls (VCF), it was found 230 times in the aligned reads (BAM), indicating that the motif is present in the mapped reads and could be biologically relevant in that specific genome region. This type of motif search can be applied in various bioinformatics analyses, such as genome annotation, regulatory element detection, and motif enrichment studies.
