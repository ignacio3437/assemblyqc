process FILTER_BAM {
    tag "$sample_id_on_tag"
    label 'process_high'

    //conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/samtools_python:97109fdca4337830"

    input:
    tuple val(sample_id_on_tag), path(hic_bam_scaffolds)

    output:
    tuple val(sample_id_on_tag), path("HiC.filtered.bam")  , emit: hic_filtered_bam_scaffolds

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    filter_bam.py $hic_bam_scaffolds 1 --NM 3 --threads $task.cpus | samtools view - -b --threads $task.cpus -o HiC.filtered.bam



    """
}
