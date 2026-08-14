#$ -N rna_flow_build_index
#$ -cwd
#$ -V
#$ -pe smp 1
#$ -l mem=16G
#$ -l h_rt=72:00:00
#$ -j y
#$ -o /home/ucsagil/Scratch/projects/rnaseq-flow/logs/rna_flow_build_index_$JOB_ID.log
#$ -m be
# send mail at begin and end of job 

PRJ_DIR=/home/ucsagil/Scratch/projects/rnaseq-flow

# Build the index for your chosen aligner
export NXF_VER=26.04.1

nextflow run ${PRJ_DIR}/rnaseq-flow/main.nf \
    --build_indices \
    --aligner star \
    --genome_fasta ${PRJ_DIR}/references/human/v102/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
    --gtf ${PRJ_DIR}/references/human/v102/Homo_sapiens.GRCh38.102.gtf.gz \
    --outdir ${PRJ_DIR}/indices/human \
    -profile singularity \
    -resume
