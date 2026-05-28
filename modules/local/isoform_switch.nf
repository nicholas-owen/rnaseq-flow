process ISOFORM_SWITCH {
    label 'process_high'
    // IsoformSwitchAnalyzeR container
    container 'quay.io/biocontainers/bioconductor-isoformswitchanalyzer:2.0.0--r43hdfd78af_0'

    input:
    path samplesheet
    path transcript_fasta
    path gtf
    path salmon_results // Per-sample Salmon output directories (named after each sample)
    path r_script       // assets/isoform_switch.R, staged so it is visible inside the container

    output:
    path "isoform_switch_output", emit: results
    path "versions.yml"         , emit: versions

    script:
    """
    Rscript ${r_script} \\
        ${samplesheet} \\
        ${transcript_fasta} \\
        ${gtf} \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-isoformswitchanalyzer: \$(Rscript -e "library(IsoformSwitchAnalyzeR); cat(as.character(packageVersion('IsoformSwitchAnalyzeR')))")
    END_VERSIONS
    """

    stub:
    """
    mkdir isoform_switch_output
    touch versions.yml
    """
}
