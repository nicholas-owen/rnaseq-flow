process RMATS {
    label 'process_high'
    container 'quay.io/biocontainers/rmats:4.3.0--py39h9ee0642_2'

    input:
    path samplesheet
    tuple val(meta), path(bams), path(bais) // Collected list of all BAMs/BAIs
    path gtf
    path script

    output:
    path "rmats_output", emit: results
    path "versions.yml", emit: versions

    script:
    // The collapsed meta carries single_end (set in workflows/rnaseq.nf).
    def run_type = meta.single_end ? "single" : "paired"
    def read_len = params.read_length ?: 100
    """
    python3 ${script} \\
        ${samplesheet} \\
        ${gtf} \\
        ${run_type} \\
        ${read_len} \\
        rmats_output \\
        ${bams}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rmats: \$(rmats.py --version 2>&1 | sed 's/v//')
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    mkdir rmats_output
    touch versions.yml
    """
}
