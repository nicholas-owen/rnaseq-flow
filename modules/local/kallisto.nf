process KALLISTO_QUANT {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/kallisto:0.48.0--h1591592_0'

    input:
    tuple val(meta), path(reads)
    path index // File (usually .idx)

    output:
    tuple val(meta), path("${meta.id}"), emit: results
    tuple val(meta), path("${meta.id}/abundance.tsv"), emit: abundance
    path "versions.yml"                , emit: versions

    script:
    // args carries the strandedness flag (--fr-stranded / --rf-stranded),
    // resolved from params.strandedness in conf/modules.config.
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Kallisto single end requires -l and -s (mean/sd fragment length), we'll assume standard or allow user args
    // Integrating via simple logic, but for robustness might need better defaults if params not set.
    def single_opts = meta.single_end ? "--single -l 200 -s 20" : ""

    """
    mkdir ${prefix}
    kallisto quant \\
        -i $index \\
        -t $task.cpus \\
        -o ${prefix} \\
        $single_opts \\
        $args \\
        $reads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallisto: \$(kallisto version | sed -e "s/kallisto, version //g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}/abundance.tsv
    touch ${prefix}/run_info.json
    touch versions.yml
    """
}
