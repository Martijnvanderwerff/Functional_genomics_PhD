#!/bin/bash
#SBATCH --job-name=cellranger_sample_1
#SBATCH --output=cellranger_sample_1.out
#SBATCH --error=cellranger_sample_1.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=100gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

ml Python

cd /groups/umcg-franke-scrna/tmp03/projects/reversetime/processed/

COMMAND=/groups/umcg-franke-scrna/tmp03/software/cellranger-arc-2.0.2/cellranger-arc
REFERENCE=/groups/umcg-franke-scrna/tmp03/external_datasets/refdata-cellranger-arc-GRCh38-2020-A-2.0.0

${COMMAND} count --id 1 \
                --libraries=/groups/umcg-franke-scrna/tmp03/projects/reversetime/raw/libraries/library_atac.csv \
                --reference=${REFERENCE} \
                --localcores=$SLURM_CPUS_PER_TASK \
                --localmem=100