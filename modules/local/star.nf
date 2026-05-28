process STAR_ALIGN {
    tag "$meta.id"
    label 'process_high'
    container 'quay.io/biocontainers/star:2.7.10b--h6b7c446_1'

    input:
    tuple val(meta), path(reads)
    path index
    path gtf

    output:
    tuple val(meta), path('*Log.final.out')        , emit: log_final
    tuple val(meta), path('*Log.out')              , emit: log_out
    tuple val(meta), path('*Log.progress.out')     , emit: log_progress
    tuple val(meta), path('*d.out.bam')            , emit: bam
    tuple val(meta), path('*ReadsPerGene.out.tab') , optional:true, emit: read_counts
    tuple val(meta), path('*Chimeric.out.junction'), optional:true, emit: chimeric_junction
    path "versions.yml"                            , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // STAR cannot read a gzipped GTF, so decompress it first if needed.
    def gtf_gz  = gtf.name.endsWith('.gz')
    def gtf_use = gtf_gz ? gtf.baseName : "${gtf}"
    // For STAR-Fusion the Chimeric.out.junction FILE is required, so chimOutType
    // must include 'Junctions'. --chimOutJunctionFormat 1 adds the summary
    // comment line that STAR-Fusion expects.
    def fusion_flags = params.ctat_lib ?
        "--chimSegmentMin 12 --chimJunctionOverhangMin 12 --chimOutType Junctions --chimOutJunctionFormat 1 --chimMultimapNmax 20" : ""
    """
    ${ gtf_gz ? "gunzip -c ${gtf} > ${gtf_use}" : "" }
    STAR \\
        --genomeDir $index \\
        --readFilesIn $reads \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix ${prefix}. \\
        --outSAMtype BAM SortedByCoordinate \\
        --readFilesCommand zcat \\
        --quantMode GeneCounts \\
        --sjdbGTFfile ${gtf_use} \\
        $fusion_flags \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.Log.final.out
    touch ${prefix}.Log.out
    touch ${prefix}.Log.progress.out
    touch ${prefix}.Aligned.sortedByCoord.out.bam
    touch ${prefix}.ReadsPerGene.out.tab
    ${ params.ctat_lib ? "touch ${prefix}.Chimeric.out.junction" : "" }
    touch versions.yml
    """
}
