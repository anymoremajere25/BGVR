nextflow.enable.dsl=2

params.bam_path = '/mnt/c/Users/trian/BGVR/chapter_06/experiment_62/src/example.bam'
params.region_list = ['1:1-50000', '1:50001-100000']

process coverageIndexing {
    input:
    val region

    output:
    file "coverage_${region.replace(':', '_')}.txt"

    script:
    def region_sanitized = region.replace(':', '_')
    """
    echo "Processing region: ${region}"
    /mnt/c/Users/trian/BGVR/chapter_06/experiment_62/target/debug/coverage_tool \
        --bam ${params.bam_path} \
        --region ${region} > coverage_${region_sanitized}.txt
    """
}

workflow {
    def regionChannel = Channel.from(params.region_list)
    coverageIndexing(regionChannel)
}