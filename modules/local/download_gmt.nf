process DOWNLOAD_GMT {
    label 'process_low'
    // No fixed container image: the r-msigdbr environment is provisioned on the
    // fly by Nextflow Wave from the Conda package below. Wave builds a minimal,
    // secure container at runtime (see the wave {} block in nextflow.config),
    // which avoids depending on a specific pre-built biocontainer tag.
    //
    // conda-forge, NOT bioconda: msigdbr is a CRAN package, so it is packaged as
    // conda-forge::r-msigdbr. It has never existed in bioconda, so the previous
    // `bioconda::r-msigdbr=7.5.1` could not be solved and this process could not
    // run at all.
    //
    // Pinned exactly, because the API is version-sensitive: msigdbr 10.0.0
    // renamed the `category` argument to `collection` and the `gs_subcat` column
    // to `gs_subcollection`, both of which assets/download_gmt.R relies on. The
    // script asserts msigdbr >= 10 at run time and fails with a clear message
    // rather than producing empty gene sets.
    //
    // From msigdbr 24.1.0 the gene sets are fetched over the network on first
    // use rather than shipped inside the package, so this process needs outbound
    // HTTPS -- see the note at the top of assets/download_gmt.R.
    conda 'conda-forge::r-msigdbr=26.1.0'

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
