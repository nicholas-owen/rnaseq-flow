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
    conda 'conda-forge::r-base=4.5.3 bioconda::bioconductor-deseq2=1.50.2 bioconda::bioconductor-apeglm=1.32.0 conda-forge::r-pheatmap=1.0.13'

    input:
    path samplesheet
    val  design
    path gene_info
    path counts_files
    path r_script
    val  ref_level
    val  padj_cutoff   // what counts as significant; see nextflow.config
    val  lfc_cutoff

    output:
    path "deseq2_output", emit: results
    path "versions.yml" , emit: versions

    script:
    // task.ext.args carries '--matrix' when the run started from a pre-computed
    // count matrix (--counts), and is empty otherwise. It is appended with its
    // own trailing space (not a separate ${} slot) so that when the flag is
    // absent the rendered command is byte-identical to the pre-feature version.
    // That keeps the task hash stable, so a normal run still resumes (-resume)
    // across this change. deseq2.R strips the flag out wherever it appears.
    def args = task.ext.args ? "${task.ext.args} " : ''
    // The baseline is exported here, from a process input, rather than set once
    // in the env scope of nextflow.config. Nextflow does not hash that scope,
    // so the old arrangement let a run with a different --reference_level
    // resume from cache and republish the previous baseline unchanged. See the
    // note where the env block used to be in nextflow.config. deseq2.R reads
    // the same variable either way, so nothing changed in the R.
    //
    // This does change the rendered command, so caches from before it are
    // invalidated once. That is the point: a different baseline is a different
    // result and should not resume.
    """
    export RNASEQ_FLOW_REFERENCE_LEVEL='${ref_level}'
    export RNASEQ_FLOW_PADJ_CUTOFF='${padj_cutoff}'
    export RNASEQ_FLOW_LFC_CUTOFF='${lfc_cutoff}'
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
