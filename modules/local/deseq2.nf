process DESEQ2 {
    label 'process_medium'
    // No fixed image: DESeq2 + apeglm are provisioned on the fly by Nextflow
    // Wave from the Conda packages below (see the wave {} block in
    // nextflow.config). lfcShrink(type="apeglm") needs the apeglm package,
    // which the deseq2-only biocontainer does not carry -- apeglm is a DESeq2
    // 'Suggests' dependency, so it is not pulled into that image.
    conda 'bioconda::bioconductor-deseq2=1.38.0 bioconda::bioconductor-apeglm'

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
    """
    Rscript ${r_script} ${samplesheet} '${design}' ${gene_info} ${counts_files}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
        apeglm: \$(Rscript -e "cat(as.character(packageVersion('apeglm')))")
    END_VERSIONS
    """

    stub:
    """
    mkdir deseq2_output
    touch versions.yml
    """
}
