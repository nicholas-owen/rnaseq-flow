process GPROFILER {
    label 'process_medium'
    // No fixed image: provisioned on the fly by Nextflow Wave from the Conda
    // packages below, as DESEQ2 and QUARTO_REPORT are.
    //
    // There is no biocontainer for gprofiler2 -- the quay.io/biocontainers
    // repository does not exist (zero tags), so the image previously named here
    // could never be pulled and this process could never run. gprofiler2 is a
    // CRAN package, hence conda-forge.
    conda 'conda-forge::r-gprofiler2=0.2.3 conda-forge::r-ggplot2'

    input:
    val organism
    path deseq2_dir
    path r_script

    output:
    path "gprofiler_output", emit: results
    path "versions.yml"    , emit: versions

    script:
    """
    Rscript ${r_script} ${organism} ${deseq2_dir}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gprofiler2: \$(Rscript -e "library(gprofiler2); cat(as.character(packageVersion('gprofiler2')))")
    END_VERSIONS
    """

    stub:
    """
    mkdir gprofiler_output
    touch versions.yml
    """
}
