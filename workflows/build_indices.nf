include { STAR_GENOME_GENERATE } from '../modules/local/star_genome_generate'
include { HISAT2_BUILD         } from '../modules/local/hisat2_build'
include { SALMON_INDEX         } from '../modules/local/salmon_index'
include { KALLISTO_INDEX       } from '../modules/local/kallisto_index'
include { GFFREAD              } from '../modules/local/gffread'

workflow BUILD_INDICES {
    take:
    aligner

    main:
    // This workflow assumes parameters provided globally via params.
    // e.g. params.genome_fasta, params.gtf, params.transcript_fasta
    
    // We check which aligner is requested (or all?)
    // User request: "create aligner specific indexes... for STAR, HISAT2, kallisto and others"
    // Usually we only build what is needed. 
    // We can use 'params.aligner' to decide, or a specific 'params.build_all_indices'.
    // Let's support individual aligner via 'params.aligner' (for building one) OR 'all'.
    
    // Pre-checks
    // Ensure Fastas are available
    // For STAR/HISAT2: params.genome_fasta, params.gtf
    // For Salmon/Kallisto: params.transcript_fasta (and genome_fasta for decoys in Salmon)
    
    // Logic:
    // If aligner == 'star' or 'all' -> Build STAR
    // If aligner == 'hisat2' or 'all' -> Build HISAT2
    // If aligner == 'salmon' or 'all' -> Build Salmon
    // If aligner == 'kallisto' or 'all' -> Build Kallisto
    
    def run_star = (aligner == 'star' || aligner == 'all')
    def run_hisat2 = (aligner == 'hisat2' || aligner == 'all')
    def run_salmon = (aligner == 'salmon' || aligner == 'all')
    def run_kallisto = (aligner == 'kallisto' || aligner == 'all')
    
    if (run_star) {
        if (!params.genome_fasta || !params.gtf) error "STAR index build requires --genome_fasta and --gtf"
        STAR_GENOME_GENERATE( file(params.genome_fasta), file(params.gtf) )
    }
    
    if (run_hisat2) {
        // GTF is required: HISAT2_BUILD extracts splice sites and exons from it
        // to build a splice-aware index.
        if (!params.genome_fasta || !params.gtf)
            error "HISAT2 index build requires --genome_fasta and --gtf"
        HISAT2_BUILD( file(params.genome_fasta), file(params.gtf) )
    }
    
    // ---- transcriptome for the pseudo-aligners ------------------------------
    // Cut from the genome with gffread rather than downloaded, so the
    // transcriptome and the annotation describe the same transcripts by
    // construction. Ensembl's own files do not, in both directions.
    //
    // cdna.all.fa.gz carries transcripts on patch and haplotype scaffolds the
    // primary-assembly GTF omits: 14212 on GRCh38.116. Every transcript-level
    // analysis had to reconcile that independently, and it is a quantification
    // problem as well as a bookkeeping one, because those haplotype copies take
    // reads from their primary counterparts.
    //
    // It also omits far more than it adds. cdna.all is cDNA only, so it holds
    // 465769 transcripts against the GTF's 646577: the missing ~198000 are
    // almost all lncRNA, which Ensembl ships separately in ncrna.fa. Most
    // lncRNAs are Pol II products, capped and poly-adenylated, so a poly-A
    // library captures them, and leaving them out of the index pushes their
    // reads onto whatever mRNA they share sequence with. Taking the whole GTF
    // is what nf-core/rnaseq does and is the reason to prefer it.
    //
    // --transcript_fasta still overrides, for a curated or non-Ensembl
    // transcriptome, but it warns: nothing can check that a supplied FASTA
    // matches the GTF, and a mismatch is silent until a downstream tool
    // refuses the pairing.
    def tx_fasta = null
    if (run_salmon || run_kallisto) {
        if (params.transcript_fasta) {
            log.warn "--transcript_fasta supplied, so the transcriptome is not derived from " +
                     "--gtf. Nothing verifies that the two describe the same transcripts, and " +
                     "a mismatch surfaces only when a transcript-level analysis rejects it. " +
                     "Omit --transcript_fasta to cut the transcriptome from the genome instead."
            tx_fasta = Channel.value(file(params.transcript_fasta))
        }
        else {
            if (!params.genome_fasta || !params.gtf)
                error "Building a ${aligner} index needs --genome_fasta and --gtf, from which " +
                      "the transcriptome is derived. Supply --transcript_fasta to use one you " +
                      "already have instead."
            GFFREAD( file(params.genome_fasta), file(params.gtf) )
            tx_fasta = GFFREAD.out.transcript_fasta
        }
    }

    if (run_salmon) {
        // The genome supplies the decoys, so Salmon needs it either way.
        if (!params.genome_fasta) error "Salmon index build with decoys requires --genome_fasta"
        SALMON_INDEX( file(params.genome_fasta), tx_fasta )
    }

    if (run_kallisto) {
        KALLISTO_INDEX( tx_fasta )
    }
}
