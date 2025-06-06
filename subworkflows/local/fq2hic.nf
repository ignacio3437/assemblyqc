include { FASTQ_FASTQC_UMITOOLS_FASTP   } from '../nf-core/fastq_fastqc_umitools_fastp/main'
include { FASTQ_BWA_MEM_SAMBLASTER      } from '../gallvp/fastq_bwa_mem_samblaster/main'
include { SEQKIT_SORT                   } from '../../modules/nf-core/seqkit/sort/main'
include { HICQC                         } from '../../modules/gallvp/hicqc'

include { BAM_FASTA_YAHS_JUICER_PRE_JUICER_TOOLS_PRE } from '../gallvp/bam_fasta_yahs_juicer_pre_juicer_tools_pre/main'

include { HIC2HTML                      } from '../../modules/local/hic2html'

workflow FQ2HIC {
    take:
    reads                           // [ val(meta), [ fq ] ]
    ref                             // [ val(meta2), fa ]
    hic_skip_fastp                  // val: true|false
    hic_skip_fastqc                 // val: true|false
    hic_alphanumeric_sort           // val: true|false

    main:
    ch_versions                     = Channel.empty()

    // SUBWORKFLOW: FASTQ_FASTQC_UMITOOLS_FASTP
    FASTQ_FASTQC_UMITOOLS_FASTP(
        reads,
        hic_skip_fastqc,
        false,                      // with_umi
        true,                       // skip_umi_extract
        0,                          // umi_discard_read
        hic_skip_fastp,
        [],                         // adapter_fasta
        true,                       // save_trimmed_fail
        false,                      // save_merged
        1                           // min_trimmed_reads
    )

    ch_fastp_log                    = FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_log
    ch_trim_reads                   = FASTQ_FASTQC_UMITOOLS_FASTP.out.reads
    ch_versions                     = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.versions)

    // MODULE: SEQKIT_SORT
    SEQKIT_SORT ( hic_alphanumeric_sort ? ref : Channel.empty() )

    ch_sorted_ref                   = hic_alphanumeric_sort
                                    ? SEQKIT_SORT.out.fastx
                                    : ref
    ch_versions                     = ch_versions.mix(SEQKIT_SORT.out.versions)

    // SUBWORKFLOW: FASTQ_BWA_MEM_SAMBLASTER
    val_sort_bam = true
    FASTQ_BWA_MEM_SAMBLASTER(
        ch_trim_reads,
        ch_sorted_ref.map { meta2, fa -> [ meta2, fa, [] ] },
        val_sort_bam
    )

    ch_bam                          = FASTQ_BWA_MEM_SAMBLASTER.out.bam
    ch_versions                     = ch_versions.mix(FASTQ_BWA_MEM_SAMBLASTER.out.versions)

    // MODULE: HICQC
    ch_bam_and_ref                  = ch_bam
                                    | map { meta, bam -> [ meta.ref_id, meta, bam ] }
                                    | join(
                                        ch_sorted_ref.map { meta2, fa -> [ meta2.id, fa ] }
                                    )
                                    | map { _ref_id, meta, bam, fa ->
                                        [ [ id: "${meta.id}.on.${meta.ref_id}", ref_id: meta.ref_id ], bam, fa ]
                                    }

    HICQC ( ch_bam_and_ref.map { meta3, bam, _fa -> [ meta3, bam ] } )

    ch_hicqc_pdf                    = HICQC.out.pdf
    ch_versions                     = ch_versions.mix(HICQC.out.versions)

    // SUBWORKFLOW: BAM_FASTA_YAHS_JUICER_PRE_JUICER_TOOLS_PRE
    BAM_FASTA_YAHS_JUICER_PRE_JUICER_TOOLS_PRE (
        ch_bam_and_ref.map { meta3, bam, _fa -> [ meta3, bam ] },
        ch_sorted_ref
    )

    ch_hic                          = BAM_FASTA_YAHS_JUICER_PRE_JUICER_TOOLS_PRE.out.hic
    ch_versions                     = ch_versions.mix(BAM_FASTA_YAHS_JUICER_PRE_JUICER_TOOLS_PRE.out.versions)

    // MODULE: HIC2HTML
    HIC2HTML ( ch_hic )

    ch_versions                     = ch_versions.mix(HIC2HTML.out.versions.first())

    emit:
    fastp_log                       = ch_fastp_log
    hicqc_pdf                       = ch_hicqc_pdf
    hic                             = ch_hic
    html                            = HIC2HTML.out.html
    versions                        = ch_versions
}
