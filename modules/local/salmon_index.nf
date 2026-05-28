process SALMON_INDEX {
    label 'process_medium'
    container 'quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0'

    input:
    path genome_fasta
    path transcript_fasta

    output:
    path "salmon_index", emit: index
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    // Decoy-aware Salmon index: the genome sequence names become decoys.
    // zcat/cat picks the right reader so this works for plain or gzipped FASTA.
    def reader = genome_fasta.name.endsWith('.gz') ? 'zcat' : 'cat'
    """
    mkdir salmon_index

    # Build the decoy list from the genome FASTA headers.
    ${reader} ${genome_fasta} | grep "^>" | cut -d " " -f 1 | sed 's/>//g' > decoys.txt

    # gentrome = transcriptome followed by genome. Concatenating two gzip
    # streams yields a valid gzip stream, so this works when both inputs are
    # gzipped (Ensembl downloads) or both are plain text.
    cat ${transcript_fasta} ${genome_fasta} > gentrome.fa.gz

    salmon index \\
        -t gentrome.fa.gz \\
        -d decoys.txt \\
        -p $task.cpus \\
        -i salmon_index \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(salmon --version | sed -e "s/salmon //g")
    END_VERSIONS
    """

    stub:
    """
    mkdir salmon_index
    touch versions.yml
    """
}
