process RMATS {
    label 'process_high'
    container 'quay.io/biocontainers/rmats:4.3.0--py39hbadf43b_5'

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
    // rMATS cannot read a compressed GTF -- it fails with "unable to parse the
    // gtf" and a UnicodeDecodeError on the gzip magic byte. Decompress first,
    // as STAR, featureCounts, GTF2BED and HISAT2_BUILD already do. The
    // reference-download workflow always produces .gz, so without this rMATS
    // fails on any normally-obtained annotation.
    def gtf_gz  = gtf.name.endsWith('.gz')
    def gtf_use = gtf_gz ? gtf.baseName : "${gtf}"
    """
    ${ gtf_gz ? "gunzip -c ${gtf} > ${gtf_use}" : "" }
    python3 ${script} \\
        ${samplesheet} \\
        ${gtf_use} \\
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
