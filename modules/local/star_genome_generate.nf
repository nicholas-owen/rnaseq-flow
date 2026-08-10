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
    // --genomeSAindexNbases must be scaled to the genome, and STAR's default of
    // 14 is sized for mammals. On a small genome STAR warns that 14 "may cause
    // seg-fault at the mapping step", so it is computed below from the actual
    // sequence length as STAR documents: min(14, log2(genomeLength)/2 - 1).
    // Yeast (12 Mb) gives 11, Drosophila 12, human 14 (i.e. unchanged).
    // An explicit --genomeSAindexNbases in ext.args suppresses the automatic
    // value rather than being passed twice.
    def sa_auto = args.contains('--genomeSAindexNbases')
        ? 'SA_FLAG=""   # --genomeSAindexNbases supplied via ext.args'
        : 'SA_FLAG="--genomeSAindexNbases ${SA_NBASES}"'
    """
    ${ fasta_gz ? "gunzip -c ${fasta} > ${fasta_use}" : "" }
    ${ gtf_gz   ? "gunzip -c ${gtf} > ${gtf_use}"     : "" }
    mkdir star_index

    GENOME_LEN=\$(grep -v '^>' ${fasta_use} | tr -d '\\n' | wc -c)
    SA_NBASES=\$(awk -v L="\$GENOME_LEN" 'BEGIN {
        n = int(log(L) / log(2) / 2 - 1)
        if (n > 14) n = 14
        if (n < 5)  n = 5
        print n
    }')
    ${sa_auto}
    echo "STAR_GENOME_GENERATE: genome length \$GENOME_LEN bp -> \$SA_FLAG"

    STAR \\
        --runMode genomeGenerate \\
        --genomeDir star_index/ \\
        --genomeFastaFiles ${fasta_use} \\
        --sjdbGTFfile ${gtf_use} \\
        \$SA_FLAG \\
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
