nextflow.enable.dsl=2

workflow HIC_PREPROCESS {
    take:
        paired_reads
    
    main:
        if (!params.hic.skip) {
            FASTP(paired_reads)
            | set { ch_cleaned_paired_reads }

            paired_reads
            .join(ch_cleaned_paired_reads, remainder: true)
            | FAST_QC
        } else {
            ch_cleaned_paired_reads = Channel.of([])
        }
    
    emit:
        cleaned_paired_reads = ch_cleaned_paired_reads
}

process FASTP {
    tag "$sample_id"
    label "process_medium"
    
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.2--h5f740d0_3':
        'quay.io/biocontainers/fastp:0.23.2--h5f740d0_3' }"

    input: 
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path('*.fastp.fastq.gz')

    script:
        """ 
        fastp \
        -i ${reads[0]} \
        -o "\$(basename ${reads[0]} .fastq.gz).fastp.fastq.gz" \
        -I ${reads[1]} \
        -O "\$(basename ${reads[1]} .fastq.gz).fastp.fastq.gz" \
        --qualified_quality_phred 20 \
        --length_required 50 \
        --thread ${task.cpus}
        """
}


process FAST_QC {
    tag "$sample_id"
    label "process_medium"
    
    publishDir "${params.outdir}/hic/fastqc", mode:'copy'
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.11.9--hdfd78af_1':
        'quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1' }"

    input:
        tuple val(sample_id), path(raw_reads), path(clean_reads)

    output:
        path '*.html'
        path '*.zip'

    script:
        """
        fastqc ${raw_reads} ${clean_reads} \
        -t ${task.cpus} \
        --nogroup
        """
}