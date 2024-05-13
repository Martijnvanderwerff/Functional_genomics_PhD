#!/bin/bash
#SBATCH --job-name=pycistopic_LDA_5000_mono
#SBATCH --output=pycistopic_LDA_5000_mono.out
#SBATCH --error=pycistopic_LDA_5000_mono.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=100gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

eval "$(conda shell.bash hook)"
conda activate scenicplus
ml Java

python3 pycistopic.py > progess.txt
