process QUARTO_REPORT {
    label 'process_single'
    // No fixed container image: the rendering environment is provisioned on the
    // fly by Nextflow Wave from the Conda packages below (see the wave {} block
    // in nextflow.config). The report needs the Quarto CLI plus R with
    // ggplot2/dplyr/tidyr/jsonlite and the interactive plotly + DT packages;
    // rocker/verse does not carry plotly/DT, so a Conda-built image is used.
    //
    // Quarto is pinned for reproducibility. The R packages are deliberately
    // left unpinned: they resolve correctly (plotly, DT and the rest are all
    // present and working) and over-constraining the environment only risks
    // solver conflicts. Pinning the whole environment belongs with publishing
    // a frozen image, which is what nf-core does: build once, push to a
    // registry, and reference it by digest rather than rebuilding per run.
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
    # The Quarto CLI is a shell wrapper that execs a bundled Deno runtime and
    # calls out to pandoc, dart-sass, esbuild and typst. It expects the upstream
    # tarball layout (a tools/x86_64/ tree beside the binary), which the Conda
    # package does not use -- it installs those tools flat into bin/ and wires
    # the paths up in an activation script instead. Nextflow runs containers
    # without activating the Conda environment, so none of that is applied and
    # the wrapper dies looking for a bundled runtime that was never there:
    #   /opt/conda/bin/quarto: line 196: .../tools/x86_64/deno: No such file
    # Sourcing the package's own activation script sets QUARTO_DENO,
    # QUARTO_PANDOC, QUARTO_DART_SASS, QUARTO_SHARE_PATH and the rest to the
    # Conda locations. Only quarto.sh is sourced: the compiler activation
    # scripts alongside it reference CONDA_PREFIX unguarded, which trips the
    # `set -u` this pipeline runs its shells with.
    _quarto_activate="\$(dirname "\$(dirname "\$(command -v quarto)")")/etc/conda/activate.d/quarto.sh"
    if [ -r "\$_quarto_activate" ]; then
        . "\$_quarto_activate"
    fi

    # Keep Deno's and Quarto's caches inside the task directory: the container
    # runs as the invoking user, whose home is not writable in the image.
    export HOME="\$PWD"
    export XDG_CACHE_HOME="\$PWD/.cache"

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
