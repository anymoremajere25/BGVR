params.bam_list = "bams.txt"
params.vcf_file = "cohort.vcf"
params.rust_bin = "/mnt/c/Users/trian/BGVR/chapter_07/experiment_71/target/debug"

process rustCoverageRunner {
    input:
    path bam_list_file
    path vcf_file
    
    output:
    path "coverage_output.txt"

    script:
    """
    BAM_FILES=\$(cat ${bam_list_file} | tr '\\n' ',' | sed 's/,\$//')
    
    ${params.rust_bin}/rust_noodles_coverage \\
        --vcf-file ${vcf_file} \\
        --bam-files "\$BAM_FILES" > coverage_output.txt
    """
}

workflow {
    bamListFile = file(params.bam_list)
    vcfFile = file(params.vcf_file)
    
    rustCoverageRunner(bamListFile, vcfFile)
}