process FEATURECOUNTS {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/subread:2.0.1--hed695b0_0'

    input:
    tuple val(meta), path(bams), path(annotation), val(strandedness)

    output:
    tuple val(meta), path("*featureCounts.txt")        , emit: counts
    tuple val(meta), path("*featureCounts.txt.summary"), emit: summary
    path "versions.yml"                                , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired = meta.single_end ? "" : "-p"

    // featureCounts strand flag: 0 = unstranded, 1 = forward, 2 = reverse.
    // 'strandedness' is the per-sample value resolved upstream (either the
    // RSeQC inference, when --strandedness auto, or the --strandedness value).
    def strand_flag = "-s 0"
    if (strandedness == 'forward') {
        strand_flag = "-s 1"
    } else if (strandedness == 'reverse') {
        strand_flag = "-s 2"
    }

    // Decompress the annotation if it is gzipped.
    def gtf_gz  = annotation.name.endsWith('.gz')
    def gtf_use = gtf_gz ? annotation.baseName : "${annotation}"
    """
    ${ gtf_gz ? "gunzip -c ${annotation} > ${gtf_use}" : "" }
    featureCounts \\
        $args \\
        $paired \\
        $strand_flag \\
        -a ${gtf_use} \\
        -o ${prefix}.featureCounts.txt \\
        -T $task.cpus \\
        ${bams}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: \$( echo \$(featureCounts -v 2>&1) | sed -e "s/featureCounts v//g" )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.featureCounts.txt
    touch ${prefix}.featureCounts.txt.summary
    touch versions.yml
    """
}
