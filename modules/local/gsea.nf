process GSEA {
    label 'process_medium'
    container 'quay.io/biocontainers/mulled-v2-ad9dd5f398966bf899ae05f8e7c54d0fb10c553f:10526027a00508e82d3e1cf2b82142279c93ea40-0'
    // Contains fgsea, dplyr, ggplot2, tibble... using a mulled container for R 4+ and bioconductor-fgsea

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
