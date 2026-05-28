process KALLISTO_INDEX {
    label 'process_medium'
    container 'quay.io/biocontainers/kallisto:0.48.0--h159158b_2'

    input:
    path transcript_fasta

    output:
    path "kallisto_index"  , emit: index
    path "versions.yml"    , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    kallisto index \\
        -i kallisto_index \\
        $transcript_fasta \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallisto: \$(kallisto version | sed -e "s/kallisto, version //g")
    END_VERSIONS
    """

    stub:
    """
    touch kallisto_index
    touch versions.yml
    """
}
