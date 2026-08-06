process GTF2BED {
    label 'process_low'
    // Not python:*-slim: those images omit procps, and Nextflow needs `ps` to
    // collect task metrics -- without it the task fails with "Command 'ps'
    // required by nextflow to collect task metrics cannot be found". The
    // biocontainers image ships it and matches the registry used elsewhere here.
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    path gtf
    path script

    output:
    path "*.bed", emit: bed
    path "versions.yml", emit: versions

    script:
    """
    python3 ${script} ${gtf} > ${gtf.baseName}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${gtf.baseName}.bed
    touch versions.yml
    """
}
