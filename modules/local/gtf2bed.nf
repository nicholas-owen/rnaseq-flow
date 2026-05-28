process GTF2BED {
    label 'process_low'
    container 'python:3.9-slim' // Helper container

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
