#!/bin/bash
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=1:00:00
#SBATCH --output=test_export.%j.out
#SBATCH --error=test_export.%j.err
#SBATCH --partition=carter-compute

echo "Job started at: $(date)"

#plink2 --pfile /carter/users/nopopko/projects/eagles/GTEx_plinkqc/data/all_hg38 \
#       --set-all-var-ids @:#:$r:$a \
#       --make-pgen \
#       --out /carter/users/nopopko/projects/eagles/GTEx_plinkqc/data/all_hg38.recode

plink2 --pfile /carter/users/nopopko/projects/eagles/GTEx_plinkqc/data/all_hg38.recode \
       --chr 1-22,X,Y \
       --make-pgen \
       --out /carter/users/nopopko/projects/eagles/GTEx_plinkqc/data/all_hg38.recode.stdchr

echo "Job ended at: $(date)"
