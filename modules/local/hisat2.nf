process HISAT2_ALIGN {
    tag "$meta.id"
    label 'process_high'
    // Mulled container bundling HISAT2 2.2.1 + samtools 1.16.1 so the BAM can be
    // sorted in the same step (the plain hisat2 container has no samtools).
    container 'quay.io/biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0'

    input:
    tuple val(meta), path(reads)
    path index // Directory containing the HISAT2 index files (*.ht2)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.log"), emit: summary
    path "versions.yml"           , emit: versions

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    // Single-end uses -U, paired-end uses -1/-2.
    def reads_arg = meta.single_end ? "-U ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    # Derive the index basename from the staged index files (genome.1.ht2 -> genome)
    INDEX=`find -L ./ -name "*.1.ht2" | sed 's/\\.1.ht2\$//'`
    [ -z "\$INDEX" ] && INDEX=`find -L ./ -name "*.1.ht2l" | sed 's/\\.1.ht2l\$//'`
    if [ -z "\$INDEX" ]; then
        echo "ERROR: HISAT2 index files (*.1.ht2) not found in the staged index directory" >&2
        exit 1
    fi

    hisat2 \\
        -x \$INDEX \\
        $reads_arg \\
        -p $task.cpus \\
        --new-summary \\
        --summary-file ${prefix}.hisat2.summary.log \\
        $args \\
        | samtools sort -@ $task.cpus -O bam -o ${prefix}.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2: \$(hisat2 --version | grep -o 'version [^ ]*' | sed 's/version //')
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.hisat2.summary.log
    touch versions.yml
    """
}
