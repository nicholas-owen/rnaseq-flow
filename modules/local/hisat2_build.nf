process HISAT2_BUILD {
    label 'process_high'
    container 'quay.io/biocontainers/hisat2:2.2.1--h1b792b2_3'

    input:
    path fasta
    path gtf      // required: used to build a splice-aware index

    output:
    path "hisat2_index", emit: index
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    // hisat2-build and the extract scripts cannot read gzipped files: decompress first.
    def fasta_gz  = fasta.name.endsWith('.gz')
    def gtf_gz    = gtf.name.endsWith('.gz')
    def fasta_use = fasta_gz ? fasta.baseName : "${fasta}"
    def gtf_use   = gtf_gz   ? gtf.baseName   : "${gtf}"
    """
    mkdir hisat2_index
    ${ fasta_gz ? "gunzip -c ${fasta} > ${fasta_use}" : "" }
    ${ gtf_gz   ? "gunzip -c ${gtf} > ${gtf_use}"     : "" }

    # Step 1: extract splice sites from the GTF
    hisat2_extract_splice_sites.py ${gtf_use} > genome.ss

    # Step 2: extract exons from the GTF
    hisat2_extract_exons.py ${gtf_use} > genome.exon

    # Step 3: build the splice-aware index using the extracted annotation
    hisat2-build \\
        -p $task.cpus \\
        --ss genome.ss \\
        --exon genome.exon \\
        $args \\
        ${fasta_use} \\
        hisat2_index/genome

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2: \$(hisat2 --version | grep -o 'version [^ ]*' | sed 's/version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir hisat2_index
    touch hisat2_index/genome.1.ht2
    touch versions.yml
    """
}
