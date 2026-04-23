#!/bin/bash
#SBATCH --mem=16G
#SBATCH -o /cellar/users/domeyer/sbatch_outs/out/%A.%x.%a.out #STDOUT
#SBATCH -e /cellar/users/domeyer/sbatch_outs/err/%A.%x.%a.err # STDERR
#SBATCH --partition=carter-compute
#SBATCH --time=336:00:00
#SBATCH --array=1-48%5

scoreScript="/cellar/users/domeyer/repos/EAGLES/workflows/bin/score_other_cohort.py"
testSets="/cellar/users/domeyer/EAGLE/cohort_test_sets"
expression_base="/cellar/users/domeyer/EAGLE/cohort_expression"

modes=(flipAllele pcregression xgb elasticnet \
        flipAlleleLDLax pcregressionLDLax xgbLDLax elasticnetLDLax \
        flipAlleleLDMed pcregressionLDMed xgbLDMed elasticnetLDMed \
        flipAlleleLDStrict pcregressionLDStrict xgbLDStrict elasticnetLDStrict)

tis="thyroid"
base_folder="thyroid_alphagenome"
cohorts=(gtex_eur tcga_eur THCA)
test_expr_tis="thyroid"
val_expr_tis="thyroid"


task=$(( $SLURM_ARRAY_TASK_ID - 1 ))
mode_idx=$(( task % 16 ))
mode=${modes[$mode_idx]}
cohort_idx=$(( task / 16 ))
cohort=${cohorts[$cohort_idx]}
    
outdir="/cellar/shared/carterlab/projects/eagle/v0.3/$base_folder/$mode"
model_folder="$outdir/models"   

    
if [[ "$mode" == elasticnet* ]] || [[ "$mode" == pcregression* ]] || [[ "$mode" == flipAllele* ]]; then
    replace_nan_arg="--replace-nan 0"
else
    replace_nan_arg=""
fi    


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