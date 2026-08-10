process GSEA {
    label 'process_medium'
    // No fixed image: provisioned on the fly by Nextflow Wave from the Conda
    // packages below, as DESEQ2 and QUARTO_REPORT are.
    //
    // This replaces a hand-written mulled-v2 hash that could not be resolved
    // from quay.io, so GSEA could never actually run. Declaring the four
    // packages assets/gsea.R loads is both verifiable and self-documenting:
    // fgsea is Bioconductor (bioconda), the rest are CRAN (conda-forge).
    conda 'bioconda::bioconductor-fgsea=1.28.0 conda-forge::r-ggplot2 conda-forge::r-dplyr conda-forge::r-tibble'

    input:
    path gmt
    path deseq2_dir
    path r_script

    output:
    path "gsea_output" , emit: results
    path "versions.yml", emit: versions

    script:
    """
    Rscript ${r_script} ${gmt} ${deseq2_dir}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgsea: \$(Rscript -e "library(fgsea); cat(as.character(packageVersion('fgsea')))")
    END_VERSIONS
    """

    stub:
    """
    mkdir gsea_output
    touch versions.yml
    """
}
