process DOWNLOAD_REFS {
    label 'process_low'
    // Use a container with python requests. 
    // quay.io/biocontainers/python-requests is not standard.
    // quay.io/biocontainers/multiqc has python/requests?
    // quay.io/biocontainers/genomepy includes python/requests.
    container 'quay.io/biocontainers/genomepy:0.16.1--pyh7cba7a3_0'

    input:
    val species
    val source
    val release   // Ensembl release to pin (e.g. 102), or 'current'
    path script

    output:
    path "references/*"    , emit: files
    path "references/*.fa.gz" , emit: fasta, optional: true // Might be unzipped? Script saves as .gz
    path "references/*.gtf.gz", emit: gtf, optional: true
    path "references/download_log.txt", emit: log
    path "versions.yml"    , emit: versions

    script:
    """
    python3 ${script} "${species}" "${source}" references "${release}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    mkdir references
    touch references/genome.fa.gz
    touch references/annotation.gtf.gz
    touch references/download_log.txt
    touch versions.yml
    """
}
