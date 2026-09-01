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
    // Authoritative sample->BAM mapping, written as a TSV for run_rmats.py.
    // meta.ids and the staged bams come from the same closure in
    // workflows/rnaseq.nf, in the same order, so pairing by index is exact.
    // The script validates the mapping against the samplesheet and fails
    // loudly on any disagreement instead of guessing from filenames.
    def bam_list = bams instanceof List ? bams : [ bams ]
    def bam_map  = [ meta.ids, bam_list ].transpose()
                       .collect { id, b -> "${id}\t${b.name}" }.join('\n')
    // rMATS cannot read a compressed GTF -- it fails with "unable to parse the
    // gtf" and a UnicodeDecodeError on the gzip magic byte. Decompress first,
    // as STAR, featureCounts, GTF2BED and HISAT2_BUILD already do. The
    // reference-download workflow always produces .gz, so without this rMATS
    // fails on any normally-obtained annotation.
    def gtf_gz  = gtf.name.endsWith('.gz')
    def gtf_use = gtf_gz ? gtf.baseName : "${gtf}"
    """
    ${ gtf_gz ? "gunzip -c ${gtf} > ${gtf_use}" : "" }
    # Quoted delimiter: sample IDs come from the samplesheet unsanitised, so a
    # '\$' or backtick in one must not be expanded by the shell. Groovy has
    # already interpolated the mapping into the block; bash must not touch it.
    cat <<'RMATS_MAP' > sample_bam.tsv
${bam_map}
RMATS_MAP

    python3 ${script} \\
        ${samplesheet} \\
        ${gtf_use} \\
        ${run_type} \\
        ${read_len} \\
        rmats_output \\
        sample_bam.tsv

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
