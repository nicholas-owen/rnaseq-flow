process ISOFORM_SWITCH {
    label 'process_high'
    // IsoformSwitchAnalyzeR container
    // 2.0.0 was never published as a biocontainer; 2.0.1 is the nearest real
    // build (a patch release of the same version).
    container 'quay.io/biocontainers/bioconductor-isoformswitchanalyzer:2.0.1--r43ha9d7317_0'

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
    // assets/isoform_switch.R hands the GTF straight to importRdata() and does
    // no gzip handling of its own, so decompress here -- the same pattern as
    // STAR, featureCounts, HISAT2_BUILD and RMATS. The reference-download
    // workflow always produces .gz.
    def gtf_gz  = gtf.name.endsWith('.gz')
    def gtf_use = gtf_gz ? gtf.baseName : "${gtf}"
    """
    ${ gtf_gz ? "gunzip -c ${gtf} > ${gtf_use}" : "" }
    Rscript ${r_script} \\
        ${samplesheet} \\
        ${transcript_fasta} \\
        ${gtf_use} \\
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
