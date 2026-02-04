process SAMTOOLS_SUBSAMPLE_SORT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam)
    val(sample_fraction)  // Fraction of reads to sample (e.g., 0.05 for 5%)

    output:
    tuple val(meta), path("*.subsampled.sorted.bam"), emit: bam
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fraction = sample_fraction ?: 0.05
    def seed = 42
    // Convert fraction to integer percentage for -s flag (0.05 -> 5, 0.10 -> 10)
    def fraction_int = (fraction * 100).toInteger()
    
    """
    # Subsample reads and sort by name in one go
    # Format for -s is SEED.FRACTION (e.g., 42.05 for 5% with seed 42)
    samtools view -s ${seed}.${fraction_int} -b ${bam} \\
        | samtools sort -n -@ ${task.cpus} ${args} -o ${prefix}.subsampled.sorted.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.subsampled.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
