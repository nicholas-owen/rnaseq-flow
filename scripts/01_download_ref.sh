#$ -N rna_flow_dwnld_ref
#$ -cwd
#$ -V
#$ -pe smp 1
#$ -l mem=4G
#$ -l h_rt=24:00:00
#$ -j y
#$ -o /home/ucsagil/Scratch/projects/rnaseq-flow/logs/rna_flow_dwnld_ref_$JOB_ID.log
#$ -m be # send mail at begin and end of job 

PRJ_DIR=/home/ucsagil/Scratch/projects/rnaseq-flow/

mkdir -p ${PRJ_DIR}/references/human

# export NXF_VER=25.10.0
nextflow run ${PRJ_DIR}/rnaseq-flow/main.nf \
    --download_refs \
    --download_species homo_sapiens \
    --download_release 102 \
    --download_gmt \
    --organism hsapiens \
    --outdir  ${PRJ_DIR}/references/human \
    -profile singularity,ucl_myriad
# genome FASTA + GTF land in references/human/v102/