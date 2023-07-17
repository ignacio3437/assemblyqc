nextflow.enable.dsl=2

process JUICER_SORT {
    tag "$sample_id_on_tag"
    label "process_high"

    input:
        tuple val(sample_id_on_tag), path(out_links_txt)

    output:
        tuple val(sample_id_on_tag), path("*sorted.links.txt"), emit: sorted_links_txt_file

    script:
        """
        sort --parallel=${task.cpus} \
        -k2,2 -k6,6 \
        $out_links_txt \
        > out.sorted.links.txt
        """
}