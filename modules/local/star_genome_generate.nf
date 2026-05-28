process STAR_GENOME_GENERATE {
    label 'process_high'
    container 'quay.io/biocontainers/star:2.7.10b--h6b7c446_1'

    input:
    path fasta
    path gtf

    output:
    path "star_index"  , emit: index
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    // STAR cannot read gzipped FASTA/GTF, so decompress them first if needed.
    def fasta_gz  = fasta.name.endsWith('.gz')
    def gtf_gz    = gtf.name.endsWith('.gz')
    def fasta_use = fasta_gz ? fasta.baseName : "${fasta}"
    def gtf_use   = gtf_gz   ? gtf.baseName   : "${gtf}"
    """
    ${ fasta_gz ? "gunzip -c ${fasta} > ${fasta_use}" : "" }
    ${ gtf_gz   ? "gunzip -c ${gtf} > ${gtf_use}"     : "" }
    mkdir star_index
    STAR \\
        --runMode genomeGenerate \\
        --genomeDir star_index/ \\
        --genomeFastaFiles ${fasta_use} \\
        --sjdbGTFfile ${gtf_use} \\
        --runThreadN $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """

    stub:
    """
    mkdir star_index
    touch star_index/SAindex
    touch versions.yml
    """
}
