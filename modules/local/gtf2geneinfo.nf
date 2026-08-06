process GTF2GENEINFO {
    label 'process_low'
    // Not python:*-slim: those images omit procps, and Nextflow needs `ps` to
    // collect task metrics (see modules/local/gtf2bed.nf).
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    path gtf
    path script

    output:
    path "gene_info.tsv", emit: gene_info
    path "versions.yml" , emit: versions

    script:
    """
    python3 ${script} ${gtf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    printf 'gene_id\\tgene_name\\tgene_biotype\\n' > gene_info.tsv
    touch versions.yml
    """
}
