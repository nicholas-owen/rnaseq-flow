process HISAT2_BUILD {
    label 'process_high'
    container 'quay.io/biocontainers/hisat2:2.2.1--h1b792b2_3'

    input:
    path fasta
    path gtf      // used to build a splice-aware index, memory permitting

    output:
    path "hisat2_index", emit: index
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''

    // Peak memory for `hisat2-build --ss --exon` is dominated by the splice graph
    // extracted from the GTF -- of the order of 200 GB for human. Rather than be
    // OOM-killed on a smaller node, fall back to a non-splice-aware index and say so.
    //
    // need_gb comes from ext.build_memory_gb (conf/base.config); task.memory is
    // what the scheduler actually granted, after resourceLimits clamped the
    // request to --max_memory.
    if (!task.memory) {
        error "[HISAT2_BUILD] No memory allocated. Configure memory for the 'process_high' label."
    }
    // ext.build_memory_gb is a plain number of gigabytes (conf/base.config), so
    // it is compared directly against what the scheduler granted -- no cast and
    // no MemoryUnit reference is needed on this side.
    def avail_gb = task.memory.toGiga()
    def need_gb  = task.ext.build_memory_gb ?: avail_gb
    def splice_aware = avail_gb >= need_gb

    // hisat2-build and the extract scripts cannot read gzipped input: decompress
    // first. The FASTA is always needed; the GTF only on the splice-aware path.
    def fasta_gz  = fasta.name.endsWith('.gz')
    def gtf_gz    = gtf.name.endsWith('.gz')
    def fasta_use = fasta_gz ? fasta.baseName : "${fasta}"
    def gtf_use   = gtf_gz   ? gtf.baseName   : "${gtf}"

    def ss = ''
    def exon = ''
    def prepare_annotation = '# Annotation not used: see the HISAT2_BUILD warning in the run log.'
    if (splice_aware) {
        log.info "[HISAT2_BUILD] ${avail_gb} GB available, ~${need_gb} GB needed: " +
                 "building a splice-aware index from ${gtf.name}."
        prepare_annotation = [
            gtf_gz ? "gunzip -c ${gtf} > ${gtf_use}" : "",
            "hisat2_extract_splice_sites.py ${gtf_use} > genome.ss",
            "hisat2_extract_exons.py ${gtf_use} > genome.exon",
        ].findAll().join('\n    ')
        ss   = '--ss genome.ss'
        exon = '--exon genome.exon'
    }
    else {
        log.warn "[HISAT2_BUILD] Only ${avail_gb} GB available but ~${need_gb} GB is needed to build a " +
                 "splice-aware index from ${gtf.name}.\n" +
                 "  Building WITHOUT --ss/--exon. HISAT2 still aligns across junctions using its own " +
                 "model, but sensitivity to novel junctions is reduced.\n" +
                 "  To get the splice-aware index, raise the ceiling with --max_memory ${need_gb}.GB " +
                 "(on a node that large).\n" +
                 "  To accept the reduced index and silence this warning, set " +
                 "--hisat2_build_memory ${avail_gb}.GB."
    }
    """
    mkdir hisat2_index
    ${ fasta_gz ? "gunzip -c ${fasta} > ${fasta_use}" : "" }
    ${prepare_annotation}

    hisat2-build \\
        -p $task.cpus \\
        $ss \\
        $exon \\
        $args \\
        ${fasta_use} \\
        hisat2_index/genome

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2: \$(hisat2 --version | grep -o 'version [^ ]*' | sed 's/version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir hisat2_index
    touch hisat2_index/genome.1.ht2
    touch versions.yml
    """
}
