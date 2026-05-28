include { STAR_GENOME_GENERATE } from '../modules/local/star_genome_generate'
include { HISAT2_BUILD         } from '../modules/local/hisat2_build'
include { SALMON_INDEX         } from '../modules/local/salmon_index'
include { KALLISTO_INDEX       } from '../modules/local/kallisto_index'

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
    
    if (run_salmon) {
        // Requires genome (for decoys) and transcript
        if (!params.transcript_fasta || !params.genome_fasta) error "Salmon index build with decoys requires --transcript_fasta AND --genome_fasta"
        SALMON_INDEX( file(params.genome_fasta), file(params.transcript_fasta) )
    }
    
    if (run_kallisto) {
        if (!params.transcript_fasta) error "Kallisto index build requires --transcript_fasta"
        KALLISTO_INDEX( file(params.transcript_fasta) )
    }
}
