nextflow.enable.dsl=2

params.sample_list = "samples.txt"
params.output_dir = "results"
params.region = "chr1:1-32"
params.mock = true  // Set to true to use mock commands instead of actual tools

workflow {
    // Define the 'samples_ch' channel inside the workflow block
    samples_ch = Channel
        .fromPath(params.sample_list)  // Load the sample list from the file
        .splitText()                   // Split by lines (text file format)
        .map { it.trim() }             // Trim whitespace and newline characters

    // Pass the samples to the alignmentOrFetch process
    bam_ch = alignmentOrFetch(samples_ch)

    // Variant calling process
    vcf_ch = variantCalling(bam_ch)

    // Merge VCF files
    mergeVcfs(vcf_ch.collect())
}

// Define the 'alignmentOrFetch' process
process alignmentOrFetch {
    input:
    val sample_id

    output:
    path "${sample_id}.bam"

    script:
    """
    # Copy BAM files from the specified path to the current working directory
    cp "/mnt/c/Users/trian/BGVR/chapter_06/experiment_64/src/${sample_id}.bam" ./
    """
}

// Define the 'variantCalling' process
process variantCalling {
    input:
    path bam

    output:
    path "variants_${bam.simpleName}.vcf"

    script:
    if (params.mock) {
        """
        # Mock variant calling - create a dummy VCF file for testing
        echo '##fileformat=VCFv4.2' > variants_${bam.simpleName}.vcf
        echo '##source=MockVariantCaller' >> variants_${bam.simpleName}.vcf
        echo '##contig=<ID=chr1,length=248956422>' >> variants_${bam.simpleName}.vcf
        echo '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${bam.simpleName}' >> variants_${bam.simpleName}.vcf
        echo 'chr1\t14653\t.\tA\tG\t100\tPASS\t.\tGT\t0/1' >> variants_${bam.simpleName}.vcf
        echo 'chr1\t14907\t.\tA\tG\t100\tPASS\t.\tGT\t0/1' >> variants_${bam.simpleName}.vcf
        echo 'chr1\t15211\t.\tG\tA\t100\tPASS\t.\tGT\t1/1' >> variants_${bam.simpleName}.vcf
        """
    } else {
        """
        # Real variant calling using GATK HaplotypeCaller
        # Make sure GATK is properly installed and in your PATH
        gatk HaplotypeCaller \
            -R /path/to/reference.fasta \
            -I ${bam} \
            -L ${params.region} \
            -O variants_${bam.simpleName}.vcf
        """
    }
}

// Define the 'mergeVcfs' process
process mergeVcfs {
    publishDir params.output_dir, mode: 'copy'

    input:
    path vcfs

    output:
    path "merged.vcf"

    script:
    if (params.mock) {
        """
        # Mock merge - concatenate the VCF files
        cat ${vcfs[0]} > header.tmp
        grep "^#" header.tmp > merged.vcf
        for vcf in ${vcfs}; do
            grep -v "^#" \$vcf >> merged.vcf || true
        done
        """
    } else {
        """
        # Real merging using bcftools
        # Make sure bcftools is properly installed and in your PATH
        bcftools merge -o merged.vcf -O v ${vcfs}
        """
    }
}
