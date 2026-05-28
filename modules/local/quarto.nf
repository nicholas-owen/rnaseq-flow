process QUARTO_REPORT {
    label 'process_single'
    // No fixed container image: the rendering environment is provisioned on the
    // fly by Nextflow Wave from the Conda packages below (see the wave {} block
    // in nextflow.config). The report needs the Quarto CLI plus R with
    // ggplot2/dplyr/tidyr/jsonlite and the interactive plotly + DT packages;
    // rocker/verse does not carry plotly/DT, so a Conda-built image is used.
    conda 'conda-forge::quarto conda-forge::r-base conda-forge::r-rmarkdown conda-forge::r-knitr conda-forge::r-ggplot2 conda-forge::r-dplyr conda-forge::r-tidyr conda-forge::r-jsonlite conda-forge::r-plotly conda-forge::r-dt'

    input:
    path "multiqc_data_dir"
    path qmd_template
    path deseq2_dir       // deseq2_output/  (or [] when DE did not run)
    path edger_dir        // edger_output/   (or [] when DE did not run)
    path gsea_dir         // gsea_output/    (or [] when GSEA did not run)
    path gprofiler_dir    // gprofiler_output/ (or [] when gProfiler did not run)

    output:
    path "analysis_report.html", emit: html
    path "versions.yml"        , emit: versions

    script:
    """
    # Render the analysis report. The .qmd discovers multiqc_data.json and the
    # deseq2/edger/gsea/gprofiler result directories staged alongside it, and
    # renders whichever sections it finds data for.
    quarto render ${qmd_template} --output analysis_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto --version)
        r-base: \$(R --version | grep "R version" | sed 's/R version //;s/ (.*//')
    END_VERSIONS
    """

    stub:
    """
    touch analysis_report.html
    touch versions.yml
    """
}
