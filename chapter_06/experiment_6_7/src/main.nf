#!/usr/bin/env nextflow

params.bam_list = 'bams.txt'
params.bcf_file = 'cohort.bcf'
params.gff_file = 'annotations.gff'
params.bam_dir = '.' // Default to current directory

// Declare the workflow
workflow {
    // Create channels from paths
    bcfChannel = Channel.fromPath(params.bcf_file)
    gffChannel = Channel.fromPath(params.gff_file)
    
    // Read the BAM list and create a channel of fully qualified file paths
    bamChannel = Channel.fromPath(params.bam_list)
                        .splitText()
                        .map { it.trim() }
                        .map { bam -> file("${params.bam_dir}/${bam}") }

    // Call processes in correct order with proper channel connections
    splitBcf(bcfChannel)
    integrateData(splitBcf.out.flatten(), bamChannel, gffChannel)
    mergeIntegrations(integrateData.out.collect())
}

process splitBcf {
    input:
    path bcf

    output:
    path "split_*.bcf"

    """
    # Index the BCF file without -t flag (use CSI index for BCF)
    bcftools index ${bcf}
    
    # Split by chromosome
    for chr in \$(bcftools index -s ${bcf} | cut -f 1); do
        bcftools view -r \$chr ${bcf} -O b -o split_\${chr}.bcf
    done
    """
}

process integrateData {
    input:
    path bcf_file
    path bam_file
    path gff_file

    output:
    path "integrated_*.json"

    """
    # Run the integration tool - use a stage flag if available,
    # otherwise ensure BAM file has an index
    if [ ! -f "${bam_file}.bai" ]; then
        # Check if samtools is available
        which samtools > /dev/null 2>&1
        if [ \$? -eq 0 ]; then
            samtools index ${bam_file}
        else
            echo "Warning: samtools not found and BAM index missing."
        fi
    fi
    
    echo "Processing BAM: ${bam_file}, BCF: ${bcf_file}, GFF: ${gff_file}"
    
    /mnt/c/Users/trian/BGVR/chapter_06/experiment_67/target/debug/rust_integrate_tool \
      --bam ${bam_file} \
      --bcf ${bcf_file} \
      --gff ${gff_file} > integrated_\$(basename ${bam_file})_\$(basename ${bcf_file}).json
    """
}

process mergeIntegrations {
    input:
    path integrated_files

    output:
    path "final_integration.json"

    """
    # Simplest approach - just concatenate the files with array brackets
    echo "[" > final_integration.json
    
    # Get count of files
    file_count=\$(echo "${integrated_files}" | wc -w)
    current=0
    
    # Process each file
    for file in ${integrated_files}; do
        current=\$((current + 1))
        cat \$file >> final_integration.json
        
        # Add comma if not the last file
        if [ \$current -lt \$file_count ]; then
            echo "," >> final_integration.json
        fi
    done
    
    echo "]" >> final_integration.json
    """
}