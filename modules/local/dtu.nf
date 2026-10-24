process DTU {
    label 'process_high'
    // tximport + DEXSeq, provisioned on the fly by Nextflow Wave from the Conda
    // packages below (no fixed container image required).
    conda 'conda-forge::r-base=4.5.3 bioconda::bioconductor-dexseq=1.56.0 bioconda::bioconductor-tximport=1.38.2'

    input:
    path samplesheet
    path quant_dirs   // per-sample Salmon/Kallisto output directories
    path gtf          // used to build the transcript-to-gene map
    val  aligner      // 'salmon' or 'kallisto'
    path r_script     // assets/dtu_dexseq.R

    output:
    path "dtu_output" , emit: results
    path "versions.yml", emit: versions

    script:
    """
    Rscript ${r_script} ${samplesheet} ${gtf} ${aligner} . ${task.cpus} ${task.memory.toGiga()}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-dexseq: \$(Rscript -e "cat(as.character(packageVersion('DEXSeq')))")
        bioconductor-tximport: \$(Rscript -e "cat(as.character(packageVersion('tximport')))")
    END_VERSIONS
    """

    stub:
    """
    mkdir dtu_output
    touch dtu_output/dtu_transcript_results.csv
    touch dtu_output/dtu_gene_qvalues.csv
    touch versions.yml
    """
}
