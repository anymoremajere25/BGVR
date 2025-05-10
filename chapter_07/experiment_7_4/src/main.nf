#!/usr/bin/env nextflow

// Enable DSL2
nextflow.enable.dsl=2

params.sample_list = 'vcf_list.txt'
params.output_vcf = 'merged.vcf'
params.threads = 4
params.format = 'bcf'
params.tool_path = "/mnt/c/Users/trian/BGVR/chapter_07/experiment_74/target/debug/rust_vcf_merge_tool"

// Define the workflow
workflow {
    // Read the list of VCF files and make them available to the process
    vcf_files = Channel.fromPath("${workflow.launchDir}/*.vcf")
                       .collect()
                       .map { files -> [file("${workflow.launchDir}/${params.sample_list}"), files] }
    
    // Execute processes
    mergedVCF = mergeVCF(vcf_files)
    generateReport(mergedVCF)
}

// Process to merge VCF files
process mergeVCF {
    cpus params.threads
    publishDir "${workflow.launchDir}/results", mode: 'copy'
    
    input:
    tuple path(vcf_list_file), path('vcf_files/*')
    
    output:
    path "merged_vcf.${params.format}"
    
    script:
    """
    # Create a local file list with the correct paths
    ls -1 vcf_files/*.vcf > local_vcf_list.txt
    
    ${params.tool_path} \\
        --vcf-list local_vcf_list.txt \\
        --out merged_vcf.${params.format} \\
        --threads ${task.cpus} \\
        --format ${params.format}
    """
}

// Process to generate report
process generateReport {
    container "python:3.9-slim"
    publishDir "${workflow.launchDir}/results", mode: 'copy'
    
    input:
    path merged_vcf
    
    output:
    path "pipeline_report.html"
    
    script:
    """
    echo "Generating report for ${merged_vcf}"
    python generate_report.py --input ${merged_vcf} || echo "<html><body><h1>Report for ${merged_vcf}</h1><p>Processing complete.</p></body></html>" > pipeline_report.html
    """
}