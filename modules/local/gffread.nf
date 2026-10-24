process GFFREAD {
    label 'process_low'
    container 'quay.io/biocontainers/gffread:0.12.7--hdcf5f25_4'

    input:
    path genome_fasta
    path gtf

    output:
    path "transcripts.fa.gz", emit: transcript_fasta
    path "versions.yml"     , emit: versions

    script:
    def args = task.ext.args ?: ''
    // gffread indexes the genome by random access and cannot read it gzipped,
    // so decompress both inputs first. Same pattern as STAR_GENOME_GENERATE.
    def fasta_gz  = genome_fasta.name.endsWith('.gz')
    def gtf_gz    = gtf.name.endsWith('.gz')
    def fasta_use = fasta_gz ? genome_fasta.baseName : "${genome_fasta}"
    def gtf_use   = gtf_gz   ? gtf.baseName         : "${gtf}"
    """
    ${ fasta_gz ? "gunzip -c ${genome_fasta} > ${fasta_use}" : "" }
    ${ gtf_gz   ? "gunzip -c ${gtf} > ${gtf_use}"            : "" }

    # -w writes the spliced transcript sequences; -g names the genome they are
    # cut from. The FASTA header is the GTF's own transcript_id, which is the
    # whole point: the transcriptome and the annotation describe the same set
    # by construction, so nothing downstream has to reconcile them.
    gffread \\
        -w transcripts.fa \\
        -g ${fasta_use} \\
        ${gtf_use} \\
        $args

    gzip -n transcripts.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """

    stub:
    """
    echo "" | gzip > transcripts.fa.gz
    touch versions.yml
    """
}
