process JUICERTOOLS_PRE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/juicertools:2.20.00--hdfd78af_0':
        'biocontainers/juicertools:2.20.00--hdfd78af_0' }"

    input:
    tuple val(meta), path(txt)
    path(sizes)

    output:
    tuple val(meta), path("*.hic"), emit: hic
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def task_mem = (task.memory.giga*0.8).intValue()
    def use_all_cpus = task.ext.use_all_cpus ? 'yes' : 'no'
	"""
    n_proc=\$(nproc 2>/dev/null || < /proc/cpuinfo grep '^process' -c)
    task_cpus=\$([ "$use_all_cpus" = "yes" ] && echo "\$n_proc" || echo "$task.cpus")

    m_proc=\$(awk '/MemTotal/ {print int(\$2 / 1024 / 1024)}' /proc/meminfo)
    avail_mem=\$([ "$use_all_cpus" = "yes" ] && echo "\$m_proc" || echo "$task_mem" )

    mkdir user_home
    export _JAVA_OPTIONS="-Djava.util.prefs.userRoot=user_prefs -Duser.home=user_home -Xms\${avail_mem}g -Xmx\${avail_mem}g"

    juicer_tools \\
        pre \\
        $args \\
        --threads $task.cpus \\
        $txt \\
        ${prefix}.hic \\
        $sizes


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        juicertools: \$(juicer_tools --version |& sed -n 's|Juicer Tools Version\\(.*\\)|\\1|p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    touch ${prefix}.hic

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        juicertools: \$(juicer_tools --version |& sed -n 's|Juicer Tools Version\\(.*\\)|\\1|p')
    END_VERSIONS
    """
}
