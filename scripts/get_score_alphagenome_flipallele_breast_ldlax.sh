#!/bin/bash
#SBATCH --mem=32G
#SBATCH -o /cellar/users/domeyer/sbatch_outs/out/%A.%x.%a.out #STDOUT
#SBATCH -e /cellar/users/domeyer/sbatch_outs/err/%A.%x.%a.err # STDERR
#SBATCH --partition=carter-compute
#SBATCH --time=336:00:00
#SBATCH --array=1-3

scoreScript="/cellar/users/domeyer/repos/EAGLES/workflows/bin/score_other_cohort.py"
testSets="/cellar/users/domeyer/EAGLE/cohort_test_sets"
expression_base="/cellar/users/domeyer/EAGLE/cohort_expression"

mode="flipAlleleLDLax"
replace_nan_arg="--replace-nan 0"

tissue="breast"
test_expr_tis="breast_mammary_tissue"
val_expr_tis="breast"
cohorts=(gtex_eur tcga_eur_female BRCA)

task=$(( $SLURM_ARRAY_TASK_ID - 1 ))
cohort=${cohorts[$task]}
    
outdir="/cellar/shared/carterlab/projects/eagle/v0.3/breast_mammary_tissue_alphagenome/$mode"
model_folder="$outdir/models"   


if [[ "$cohort" == gtex* ]]; then
    expr_arg="--expr ${expression_base}/${cohort}/${test_expr_tis}.tsv"
    pbase="$testSets"
elif [[ "$cohort" == tcga* ]]; then
    expr_arg="--expr ${expression_base}/${cohort}/${val_expr_tis}.tsv"
    pbase="$testSets"
elif [[ -f "$testSets/tcga_tumor_sets/${cohort}.pkl" ]]; then
    expr_arg=""
    pbase="$testSets/tcga_tumor_sets"
else
    echo "cohort must be tcga/gtex/geuvadis or tumor. $cohort given." >&2
    exit 1
fi


if [ -d "$model_folder" ]; then
    score_outdir="$outdir/${cohort}"
    mkdir -p $score_outdir

    conda run -n eagle python $scoreScript \
        --pgen "$pbase/$cohort.pgen" \
        --pvar "$pbase/$cohort.pvar" \
        --psam "$pbase/$cohort.psam" \
        --pindex "$pbase/$cohort.pkl" \
        --model-folder  $model_folder \
        --outdir $score_outdir     \
        $replace_nan_arg    \
        $expr_arg 
else
    echo "Directory $model_folder does not exist"
fi