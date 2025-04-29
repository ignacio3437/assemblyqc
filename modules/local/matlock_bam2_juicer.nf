process MATLOCK_BAM2_JUICER {
    tag "$sample_id_on_tag"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/59/59df5236cc790cac47380a5654f39052a1cb3f9c7868ed397e4b3205a9fb2776/data':
        'community.wave.seqera.io/library/matlock_samtools:3c30bc2808902fde' }"

    input:
    tuple val(sample_id_on_tag), path(hic_bam_scaffolds)

    output:
    tuple val(sample_id_on_tag), path("out.links.txt")  , emit: links
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def VERSION = '20181227'
    def args2 = task.ext.args2 ?: ''
    """
    samtools \\
        view \\
        $args2 \\
        -Sb \\
        $hic_bam_scaffolds \\
        > ${sample_id_on_tag}.juicer.bam

    matlock \\
        bam2 \\
        juicer \\
        ${sample_id_on_tag}.juicer.bam \\
        out.links.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        matlock: $VERSION
    END_VERSIONS
    """

    stub:
    def VERSION = '20181227'
    """
    touch out.links.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        matlock: $VERSION
    END_VERSIONS
    """
}
