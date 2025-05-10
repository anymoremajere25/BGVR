params.bam_list     = 'bams.txt'
params.region_list  = 'regions.txt'
params.outdir       = '.'

// Add a process to index the BAM files if needed
process indexBam {
    input:
    path bam_file

    output:
    tuple path(bam_file), path("*.bai"), emit: indexed_bam

    script:
    """
    samtools index ${bam_file}
    """
}

process coverageCalc {
    tag "${bam_file}:${region}"
    
    input:
    tuple path(bam_file), path(bam_index)
    val region

    output:
    path "coverage_*.txt", emit: coverage_file

    script:
    def bam_basename = bam_file.name
    def region_sanitized = region.toString().replaceAll(/[:\-]/, '_').trim()
    def output_name = "coverage_${bam_basename}_${region_sanitized}.txt"

    """
    # Run your coverage tool
    /mnt/c/Users/trian/BGVR/chapter_06/experiment_68/target/debug/rust_coverage_tool \\
      --bam ${bam_file} \\
      --region ${region} \\
      --out ${output_name}
    
    # Make sure the output file exists no matter what
    touch ${output_name}
    """
}

process mergeCoverage {
    publishDir params.outdir, mode: 'copy'
    
    input:
    path('coverage_*')

    output:
    path "merged_coverage.txt", emit: merged_file

    script:
    """
    echo "# Merged coverage results" > merged_coverage.txt
    cat coverage_* >> merged_coverage.txt
    """
}

workflow {
    // Read BAM files as paths
    bamFiles = Channel.fromPath(params.bam_list)
        .splitText()
        .map { it.trim() }
        .filter { it.length() > 0 }
        .map { file(it.trim()) }

    // Read regions as strings
    regions = Channel.fromPath(params.region_list)
        .splitText()
        .map { it.trim() }
        .filter { it.length() > 0 }

    // Index BAM files
    indexedBams = indexBam(bamFiles)
    
    // Create combinations
    combinations = indexedBams.combine(regions)

    // Run coverage calculation
    coverageResults = coverageCalc(combinations.map { tuple(it[0], it[1]) }, combinations.map { it[2] })

    // Merge results
    mergedResults = mergeCoverage(coverageResults.coverage_file.collect())
    
    // You can use this to see the final output path
    mergedResults.merged_file.view { "Final merged output: $it" }
}