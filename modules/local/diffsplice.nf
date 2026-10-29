process DIFFSPLICE {
    label 'process_high'
    // edgeR (+ tximport for the transcript level), provisioned on the fly by
    // Nextflow Wave from the Conda packages below - no fixed image required.
    conda 'conda-forge::r-base=4.5.3 bioconda::bioconductor-edger=4.8.2 bioconda::bioconductor-tximport=1.38.2'

    input:
    path samplesheet
    path count_input   // exon featureCounts files OR per-sample quant directories
    path gtf           // transcript-to-gene map (used by the transcript level)
    val  level         // 'exon' or 'transcript'
    val  aligner       // star | hisat2 | salmon | kallisto
    path r_script      // assets/diffsplice_edger.R
    val  ref_level     // baseline condition; see nextflow.config

    output:
    path "diffsplice_output", emit: results
    path "versions.yml"     , emit: versions

    script:
    // Baseline exported from a process input so it reaches the task hash. See
    // modules/local/deseq2.nf and the note in nextflow.config.
    """
    export RNASEQ_FLOW_REFERENCE_LEVEL='${ref_level}'
    Rscript ${r_script} ${samplesheet} ${level} ${gtf} ${aligner} . ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-edger: \$(Rscript -e "cat(as.character(packageVersion('edgeR')))")
    END_VERSIONS
    """

    stub:
    """
    mkdir diffsplice_output
    touch diffsplice_output/diffsplice_${level}_results.csv
    touch diffsplice_output/diffsplice_gene_results.csv
    touch versions.yml
    """
}
