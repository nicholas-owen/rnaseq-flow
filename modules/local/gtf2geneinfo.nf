process GTF2GENEINFO {
    label 'process_low'
    container 'python:3.9-slim' // Helper container

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
