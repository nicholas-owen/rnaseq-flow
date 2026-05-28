process RSEQC {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/rseqc:5.0.1--py38h24c8516_2'

    input:
    tuple val(meta), path(bam), path(bai)
    path bed           // BED12 annotation
    path strand_script // assets/parse_strandedness.py

    output:
    tuple val(meta), path("*.infer_experiment.txt") , emit: infer_experiment
    tuple val(meta), path("*.read_distribution.txt"), emit: read_distribution
    tuple val(meta), path("*.geneBodyCoverage.txt") , emit: genebody_coverage
    tuple val(meta), path("*.strandedness.txt")     , emit: strandedness
    tuple val(meta), path("*.pdf")                  , emit: pdfs
    path "versions.yml"                             , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # 1. Infer experiment (library strandedness)
    infer_experiment.py -i ${bam} -r ${bed} > ${prefix}.infer_experiment.txt

    # 2. Read distribution over genomic features
    read_distribution.py -i ${bam} -r ${bed} > ${prefix}.read_distribution.txt

    # 3. Gene body coverage (5'-3' evenness)
    geneBody_coverage.py -i ${bam} -r ${bed} -o ${prefix}

    # 4. Reduce infer_experiment to a single strandedness verdict so downstream
    #    steps (featureCounts) can consume it directly.
    python3 ${strand_script} ${prefix}.infer_experiment.txt > ${prefix}.strandedness.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(infer_experiment.py --version | sed 's/infer_experiment.py //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.infer_experiment.txt
    touch ${prefix}.read_distribution.txt
    touch ${prefix}.geneBodyCoverage.txt
    touch ${prefix}.geneBodyCoverage.curves.pdf
    echo unstranded > ${prefix}.strandedness.txt
    touch versions.yml
    """
}
