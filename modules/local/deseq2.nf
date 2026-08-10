process DESEQ2 {
    label 'process_medium'
    // No fixed image: DESeq2 + apeglm are provisioned on the fly by Nextflow
    // Wave from the Conda packages below (see the wave {} block in
    // nextflow.config). lfcShrink(type="apeglm") needs the apeglm package,
    // which the deseq2-only biocontainer does not carry -- apeglm is a DESeq2
    // 'Suggests' dependency, so it is not pulled into that image.
    //
    // r-pheatmap draws the top-variable-gene heatmap. deseq2.R guards its use
    // with requireNamespace() and skips the heatmap when it is absent, so
    // omitting it here did not fail the run -- it silently produced no
    // heatmap_top_var.png while the report still described one. It is left
    // unpinned so the solver can pick the build matching the R version that
    // DESeq2 1.38.0 pulls in.
    conda 'bioconda::bioconductor-deseq2=1.38.0 bioconda::bioconductor-apeglm conda-forge::r-pheatmap'

    input:
    path samplesheet
    val  design
    path gene_info
    path counts_files
    path r_script

    output:
    path "deseq2_output", emit: results
    path "versions.yml" , emit: versions

    script:
    // task.ext.args carries '--matrix' when the run started from a pre-computed
    // count matrix (--counts), and is empty otherwise. It is appended with its
    // own trailing space (not a separate ${} slot) so that when the flag is
    // absent the rendered command is byte-identical to the pre-feature version
    // — that keeps the task hash stable, so a normal run still resumes (-resume)
    // across this change. deseq2.R strips the flag out wherever it appears.
    def args = task.ext.args ? "${task.ext.args} " : ''
    """
    Rscript ${r_script} ${samplesheet} '${design}' ${gene_info} ${args}${counts_files}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
        apeglm: \$(Rscript -e "cat(as.character(packageVersion('apeglm')))")
        pheatmap: \$(Rscript -e "cat(as.character(packageVersion('pheatmap')))" 2>/dev/null || echo "not installed")
    END_VERSIONS
    """

    stub:
    """
    mkdir deseq2_output
    touch versions.yml
    """
}
