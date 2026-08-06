process QUARTO_REPORT {
    label 'process_single'
    // No fixed container image: the rendering environment is provisioned on the
    // fly by Nextflow Wave from the Conda packages below (see the wave {} block
    // in nextflow.config). The report needs the Quarto CLI plus R with
    // ggplot2/dplyr/tidyr/jsonlite and the interactive plotly + DT packages;
    // rocker/verse does not carry plotly/DT, so a Conda-built image is used.
    //
    // Quarto is pinned to 1.7.31. Left unpinned it resolved to 1.9.38, whose
    // Wave-built image is incomplete: the Quarto CLI is a shell wrapper that
    // execs a bundled Deno runtime, and the build ships neither
    // bin/tools/x86_64/deno (only typst-gather survives there) nor
    // share/version, so `quarto render` dies with
    //   /opt/conda/bin/quarto: line 208: .../tools/x86_64/deno: No such file
    // The R side of the same image is complete, so this is specific to how that
    // Quarto version is packaged/pruned. 1.7.31 is the version nf-core's
    // quartonotebook module ships, in an image built the same way, so it is
    // known to survive this route.
    //
    // The R packages are deliberately left unpinned: they resolved correctly
    // (plotly, DT and the rest were all present and working) and constraining
    // them alongside an older Quarto only risks solver conflicts. Pinning the
    // whole environment belongs with publishing a frozen image -- which is the
    // real fix, and what nf-core does: build once, push to a registry, and
    // reference it by digest rather than rebuilding per run.
    conda 'conda-forge::quarto=1.7.31 conda-forge::r-base conda-forge::r-rmarkdown conda-forge::r-knitr conda-forge::r-ggplot2 conda-forge::r-dplyr conda-forge::r-tidyr conda-forge::r-jsonlite conda-forge::r-plotly conda-forge::r-dt'

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
