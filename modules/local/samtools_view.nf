nextflow.enable.dsl=2

process SAMTOOLS_VIEW {
    tag "$sample_id_on_tag"
    label "process_high"
    
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1':
        'quay.io/biocontainers/samtools:1.16.1--h6899075_1' }"

    input:
        tuple val(sample_id_on_tag), path(marked_sam)

    output:
        tuple val(sample_id_on_tag), path("*_dedup.bam"), emit: dedup_bam

    script:
        """
        file_name="${marked_sam}"
        samtools view \
        --threads ${task.cpus} \
        -S \
        -b \
        -h \
        -F 2316 \
        "${marked_sam}" \
        > "\${file_name%%.*}_dedup.bam"
        """
}