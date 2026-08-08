/*
 * FastQC: Quality control for raw sequencing reads
 */
process FASTQC {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip
    path  "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    // FastQC is run on the reads under their own filenames, deliberately: the
    // published fastqc/<accession>_fastqc.html and the paths in
    // multiqc_sources.txt are the run's provenance, and for a public dataset the
    // accession is the identifier that ties results back to the archive.
    //
    // That does mean FastQC, fastp and STAR each report the same sample under a
    // different name, which MultiQC would otherwise be unable to merge. The
    // reconciliation is done in the reporting layer instead -- MULTIQC is given
    // a sample-name mapping built from the samplesheet (see modules/local/
    // multiqc.nf) -- so the display is normalised without rewriting any
    // filenames.
    """
    fastqc $args --quiet --threads $task.cpus $reads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed -e "s/FastQC v//g" )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_fastqc.html
    touch ${prefix}_fastqc.zip
    touch versions.yml
    """
}
