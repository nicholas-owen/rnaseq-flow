process ISOFORM_SWITCH {
    label 'process_high'
    // No fixed container image: the environment is provisioned on the fly by
    // Nextflow Wave from the Conda packages below, as DESEQ2, DTU, DIFFSPLICE
    // and the other R processes already are.
    //
    // This used to pin quay.io/biocontainers/bioconductor-isoformswitchanalyzer
    // 2.10.0--r45h01b2380_0. That image has no MASS: bioconda's r-base omits
    // R's recommended packages, and importRdata needs MASS to batch-correct
    // once it detects a covariate, so --isoform_switch died at step 6 of 10
    // with "MASS package required but is not installed".
    //
    // Declaring conda alongside the container would not have fixed it. Wave's
    // default strategy resolves 'container' before 'conda', so the spec would
    // have been ignored and the image used unchanged.
    //
    // The solve returns the same IsoformSwitchAnalyzeR build the biocontainer
    // carried, r45h01b2380_0, so this changes how the environment is assembled
    // rather than which software runs.
    conda 'conda-forge::r-base=4.5.3 bioconda::bioconductor-isoformswitchanalyzer=2.10.0 conda-forge::r-mass=7.3_66'

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
