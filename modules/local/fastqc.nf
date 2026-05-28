/*
 * FastQC: Quality control for raw sequencing reads
 */
process FASTQC {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip
    path  "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    fastqc $args --quiet --threads $task.cpus $reads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed -e "s/FastQC v//g" )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_fastqc.html
    touch ${prefix}_fastqc.zip
    touch versions.yml
    """
}
