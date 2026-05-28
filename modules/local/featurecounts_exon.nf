process FEATURECOUNTS_EXON {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/subread:2.0.1--hed695b0_0'

    input:
    tuple val(meta), path(bam), path(annotation), val(strandedness)

    output:
    tuple val(meta), path("*.exon.featureCounts.txt"), emit: counts
    path "versions.yml"                              , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired = meta.single_end ? "" : "-p"

    // featureCounts strand flag: 0 = unstranded, 1 = forward, 2 = reverse.
    def strand_flag = "-s 0"
    if (strandedness == 'forward') {
        strand_flag = "-s 1"
    } else if (strandedness == 'reverse') {
        strand_flag = "-s 2"
    }

    // Decompress the annotation if it is gzipped.
    def gtf_gz  = annotation.name.endsWith('.gz')
    def gtf_use = gtf_gz ? annotation.baseName : "${annotation}"
    // -f counts at the individual feature (exon) level instead of summarising
    // to genes; -O lets a read overlapping several exons count for each. The
    // result is a per-exon count matrix carrying the gene_id, the input the
    // edgeR diffSplice exon-usage test needs.
    """
    ${ gtf_gz ? "gunzip -c ${annotation} > ${gtf_use}" : "" }
    featureCounts \\
        -f -O \\
        $paired \\
        $strand_flag \\
        -t exon \\
        -g gene_id \\
        -a ${gtf_use} \\
        -o ${prefix}.exon.featureCounts.txt \\
        -T $task.cpus \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: \$( echo \$(featureCounts -v 2>&1) | sed -e "s/featureCounts v//g" )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.exon.featureCounts.txt
    touch versions.yml
    """
}
