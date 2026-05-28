process STAR_FUSION {
    tag "$meta.id"
    label 'process_high'
    container 'quay.io/biocontainers/star-fusion:1.12.0--hdfd78af_1'

    input:
    tuple val(meta), path(junctions) // Chimeric.out.junction from STAR_ALIGN
    path ctat_lib                    // CTAT genome library directory

    output:
    path "*.star-fusion.fusion_predictions.tsv"         , emit: fusions
    path "*.star-fusion.fusion_predictions.abridged.tsv", emit: fusions_abridged
    path "versions.yml"                                 , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // We run STAR-Fusion in its "kickstart" mode using the Chimeric.out.junction
    // file already produced by STAR_ALIGN ( -J ). FASTQ input is not required
    // in this mode. The CTAT library is expected to be an uncompressed directory.
    """
    STAR-Fusion \\
        --genome_lib_dir $ctat_lib \\
        -J $junctions \\
        --CPU $task.cpus \\
        --output_dir ${prefix}_star_fusion_out

    # Prefix outputs with the sample id so they are unique when collected.
    mv ${prefix}_star_fusion_out/star-fusion.fusion_predictions.tsv \\
       ${prefix}.star-fusion.fusion_predictions.tsv
    mv ${prefix}_star_fusion_out/star-fusion.fusion_predictions.abridged.tsv \\
       ${prefix}.star-fusion.fusion_predictions.abridged.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star_fusion: \$(STAR-Fusion --version 2>&1 | grep 'STAR-Fusion version' | sed 's/STAR-Fusion version: //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.star-fusion.fusion_predictions.tsv
    touch ${prefix}.star-fusion.fusion_predictions.abridged.tsv
    touch versions.yml
    """
}
