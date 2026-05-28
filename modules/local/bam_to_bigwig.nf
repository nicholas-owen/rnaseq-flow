process BAM_TO_BIGWIG {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/deeptools:3.5.2--pydaf84f49_1'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.bw"), emit: bigwig
    path "versions.yml"          , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Handle strandedness for BigWig? 
    // Usually standard tracks are unstranded. To separate strands, we need two runs with --filterRNAstrand.
    // For now, let's output a single unstranded coverage track for simplicity, unless requested.
    // We normalize using CPM for comparability.
    """
    bamCoverage \\
        --bam $bam \\
        --outFileName ${prefix}.bw \\
        --normalizeUsing CPM \\
        --numberOfProcessors $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(bamCoverage --version | sed 's/bamCoverage //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bw
    touch versions.yml
    """
}
