process DOWNLOAD_GMT {
    label 'process_low'
    // No fixed container image: the r-msigdbr environment is provisioned on the
    // fly by Nextflow Wave from the Conda package below. Wave builds a minimal,
    // secure container at runtime (see the wave {} block in nextflow.config),
    // which avoids depending on a specific pre-built biocontainer tag.
    conda 'conda-forge::r-msigdbr=7.5.1'

    input:
    val organism
    path script

    output:
    path "gmt/*.gmt"   , emit: gmt_files
    path "versions.yml", emit: versions

    script:
    // 'organism' may be a gProfiler-style id (e.g. hsapiens) or a scientific
    // name; download_gmt.R maps common ids to the scientific name msigdbr needs.
    """
    Rscript ${script} "${organism}" gmt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-msigdbr: \$(Rscript -e "library(msigdbr); cat(as.character(packageVersion('msigdbr')))")
    END_VERSIONS
    """

    stub:
    """
    mkdir gmt
    touch gmt/hallmark.gmt
    touch versions.yml
    """
}
