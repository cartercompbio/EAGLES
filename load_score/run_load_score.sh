#! /bin/bash
#SBATCH --mem=64G
#SBATCH -o /cellar/users/domeyer/sbatch_outs/out/%A.%x.%a.out #STDOUT
#SBATCH -e /cellar/users/domeyer/sbatch_outs/err/%A.%x.%a.err # STDERR
#SBATCH --partition=carter-compute


conda run -n eagle nextflow /cellar/users/domeyer/EAGLE/load_score/load_score.nf
