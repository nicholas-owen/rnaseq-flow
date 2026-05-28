process SALMON_QUANT {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/salmon:1.10.0--h7e5ed60_0'

    input:
    tuple val(meta), path(reads)
    path index // Directory

    output:
    tuple val(meta), path("${meta.id}"), emit: results
    tuple val(meta), path("${meta.id}/quant.sf"), emit: quant
    path "versions.yml"                , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    // Salmon requires determining library type or -l A
    if (meta.single_end) {
        """
        salmon quant \\
            --index $index \\
            --libType A \\
            --unmated_reads $reads \\
            --threads $task.cpus \\
            --output $prefix \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            salmon: \$(salmon --version | sed -e "s/salmon //g")
        END_VERSIONS
        """
    } else {
        """
        salmon quant \\
            --index $index \\
            --libType A \\
            --mates1 ${reads[0]} --mates2 ${reads[1]} \\
            --threads $task.cpus \\
            --output $prefix \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            salmon: \$(salmon --version | sed -e "s/salmon //g")
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}/quant.sf
    touch versions.yml
    """
}
