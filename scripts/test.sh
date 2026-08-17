#$ -N rna_flow_test
#$ -cwd
#$ -V
#$ -pe smp 2
#$ -l mem=12G
#$ -l h_rt=24:00:00
#$ -l tmpfs=10G
#$ -j y
#$ -o /home/ucsagil/Scratch/projects/rnaseq-flow/logs/rna_flow_test_$JOB_ID.log
#$ -m be # send mail at begin and end of job 

### ARGS for download ref, build index, run analysis. eventually to be command line controlled
run_download=F

# aligners
run_star=F
run_hisat2=F
run_salmon=T
run_kallisto=T

PRJ_DIR=/home/ucsagil/Scratch/projects/rnaseq-flow

export SINGULARITY_TMPDIR="$TMPDIR"
export APPTAINER_TMPDIR="$TMPDIR"

mkdir -p ${PRJ_DIR}/references/yeast

if [[ "$run_download" == "T" ]]; then
# export NXF_VER=25.10.0
    nextflow run ${PRJ_DIR}/rnaseq-flow/main.nf --download_refs \
        --download_species saccharomyces_cerevisiae --download_release 116 \
        --outdir ${PRJ_DIR}/references/yeast \
        -profile test_yeast,singularity \
        --max_cpus 2 \
        --max_memory 24GB
        # --download_gmt \

fi
    
if [[ "$run_star" == "T" ]]; then

    if [[ ! -d ${PRJ_DIR}/references/idx/yeast_star/star_index ]]; then
        nextflow run ${PRJ_DIR}/rnaseq-flow/main.nf --build_indices --aligner star \
            --genome_fasta ${PRJ_DIR}/references/yeast/v116/*.dna.toplevel.fa.gz \
            --gtf          ${PRJ_DIR}/references/yeast/v116/*.gtf.gz \
            --outdir ${PRJ_DIR}/references/idx/yeast_star \
            -profile test_yeast,singularity \
            --max_cpus 2 \
            --max_memory 24GB  \
        || { echo "STAR index build failed (exit $?) - skipping alignment"; exit 1; }
                else
            echo "STAR index present - skipping build"
        fi

    nextflow run ${PRJ_DIR}/rnaseq-flow/main.nf \
        --input      ${PRJ_DIR}/test_data/samplesheet_yeast.csv \
        --aligner    star \
        --star_index ${PRJ_DIR}/references/idx/yeast_star/star_index \
        --gtf        ${PRJ_DIR}/references/yeast/v116/*.gtf.gz \
        --outdir     ${PRJ_DIR}/results/test/star \
        -profile test_yeast,singularity \
        --max_cpus 2 \
        --max_memory 24GB \
        -resume
        # --gmt        ${PRJ_DIR}/references/yeast/gmt/c5_go_bp.gmt \

fi

if [[ "$run_hisat2" == "T" ]]; then

    if [[ ! -d ${PRJ_DIR}/references/idx/yeast_hisat2/hisat2_index ]]; then
        nextflow run ${PRJ_DIR}/rnaseq-flow/main.nf --build_indices --aligner hisat2 \
            --genome_fasta ${PRJ_DIR}/references/yeast/v116/*.dna.toplevel.fa.gz \
            --gtf          ${PRJ_DIR}/references/yeast/v116/*.gtf.gz \
            --outdir ${PRJ_DIR}/references/idx/yeast_hisat2/hisat2_index \
            -profile test_yeast,singularity \
            --max_cpus 2 \
            --max_memory 24GB  \
        || { echo "hisat index build failed (exit $?) - skipping alignment"; exit 1; }
                else
            echo "hisat index present - skipping build"
        fi

    nextflow run ${PRJ_DIR}/rnaseq-flow/main.nf \
        --input      ${PRJ_DIR}/test_data/samplesheet_yeast.csv \
        --aligner    hisat2 \
        --hisat2_index ${PRJ_DIR}/references/idx/yeast_hisat2/hisat2_index \
        --strandedness unstranded \
        --gtf        ${PRJ_DIR}/references/yeast/v116/*.gtf.gz \
        --outdir     ${PRJ_DIR}/results/test/hisat2 \
        -profile test_yeast,singularity \
        --max_cpus 2 \
        --max_memory 24GB \
        -resume
        # --gmt        ${PRJ_DIR}/references/yeast/gmt/c5_go_bp.gmt \

fi


if [[ "$run_salmon" == "T" ]]; then

    if [[ ! -d ${PRJ_DIR}/references/idx/yeast_salmon/salmon_index ]]; then
        nextflow run ${PRJ_DIR}/rnaseq-flow/main.nf --build_indices --aligner salmon \
            --genome_fasta ${PRJ_DIR}/references/yeast/v116/*.dna.toplevel.fa.gz \
            --transcript_fasta ${PRJ_DIR}/references/yeast/v116/*.cdna.all.fa.gz \
            --outdir ${PRJ_DIR}/references/idx/yeast_salmon \
            -profile test_yeast,singularity \
            --max_cpus 2 \
            --max_memory 24GB  \
        || { echo "salmon index build failed (exit $?) - skipping alignment"; exit 1; }
                else
            echo "salmon index present - skipping build"
        fi

    nextflow run ${PRJ_DIR}/rnaseq-flow/main.nf \
        --input      ${PRJ_DIR}/test_data/samplesheet_yeast.csv \
        --aligner    salmon \
        --salmon_index ${PRJ_DIR}/references/idx/yeast_salmon/salmon_index \
        --strandedness unstranded \
        --gtf        ${PRJ_DIR}/references/yeast/v116/*.gtf.gz \
        --outdir     ${PRJ_DIR}/results/test/salmon \
        -profile test_yeast,singularity \
        --max_cpus 2 \
        --max_memory 24GB \
        -resume
        # --gmt        ${PRJ_DIR}/references/yeast/gmt/c5_go_bp.gmt \

fi

if [[ "$run_kallisto" == "T" ]]; then

    if [[ ! -d ${PRJ_DIR}/references/idx/yeast_kallisto/kallisto_index ]]; then
        nextflow run ${PRJ_DIR}/rnaseq-flow/main.nf --build_indices --aligner kallisto \
            --transcript_fasta ${PRJ_DIR}/references/yeast/v116/*.cdna.all.fa.gz \
            --outdir ${PRJ_DIR}/references/idx/yeast_kallisto \
            -profile test_yeast,singularity \
            --max_cpus 2 \
            --max_memory 24GB  \
        || { echo "kallisto index build failed (exit $?) - skipping alignment"; exit 1; }
                else
            echo "kallisto index present - skipping build"
        fi

    nextflow run ${PRJ_DIR}/rnaseq-flow/main.nf \
        --input      ${PRJ_DIR}/test_data/samplesheet_yeast.csv \
        --aligner    kallisto \
        --kallisto_index ${PRJ_DIR}/references/idx/yeast_kallisto/kallisto_index \
        --strandedness unstranded \
        --gtf        ${PRJ_DIR}/references/yeast/v116/*.gtf.gz \
        --outdir     ${PRJ_DIR}/results/test/kallisto \
        -profile test_yeast,singularity \
        --max_cpus 2 \
        --max_memory 24GB \
        -resume
        # --gmt        ${PRJ_DIR}/references/yeast/gmt/c5_go_bp.gmt \

fi