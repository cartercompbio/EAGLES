#!/bin/bash
#SBATCH --mem=16G
#SBATCH -o /cellar/users/domeyer/sbatch_outs/out/%A.%x.%a.out #STDOUT
#SBATCH -e /cellar/users/domeyer/sbatch_outs/err/%A.%x.%a.err # STDERR
#SBATCH --partition=carter-compute
#SBATCH --time=336:00:00
#SBATCH --array=1-32%16

scoreScript="/cellar/users/domeyer/repos/EAGLES/workflows/bin/score_other_cohort.py"
testSets="/cellar/users/domeyer/EAGLE/cohort_test_sets"
expression_base="/cellar/users/domeyer/EAGLE/cohort_expression"

modes=(flipAllele pcregression elasticnet xgb \
        flipAlleleLDLax pcregressionLDLax elasticnetLDLax xgbLDLax \
        flipAlleleLDMed pcregressionLDMed elasticnetLDMed xgbLDMed \
        flipAlleleLDStrict pcregressionLDStrict elasticnetLDStrict xgbLDStrict)

tissues=(colon_transverse colon_sigmoid)
cohorts=(gtex_eur_transverse gtex_eur_sigmoid)

task=$(( $SLURM_ARRAY_TASK_ID - 1 ))
mode_idx=$(( task % 16 ))
mode=${modes[$mode_idx]}
tis_idx=$(( task / 16 ))
tis=${tissues[$tis_idx]}
cohort=${cohorts[$tis_idx]}

expr_arg="--expr ${expression_base}/gtex_eur/${tis}.tsv"
pbase="$testSets"
    
outdir="/cellar/shared/carterlab/projects/eagle/v0.3/colon_alphagenome/$mode"
model_folder="$outdir/models"   

    
if [[ "$mode" == elasticnet* ]] || [[ "$mode" == pcregression* ]] || [[ "$mode" == flipAllele* ]]; then
    replace_nan_arg="--replace-nan 0"
else
    replace_nan_arg=""
fi    

if [ -d "$model_folder" ]; then
    score_outdir="$outdir/${cohort}"
    mkdir -p $score_outdir

    conda run -n eagle python $scoreScript \
        --pgen "$pbase/gtex_eur.pgen" \
        --pvar "$pbase/gtex_eur.pvar" \
        --psam "$pbase/gtex_eur.psam" \
        --pindex "$pbase/gtex_eur.pkl" \
        --model-folder  $model_folder \
        --outdir $score_outdir     \
        $replace_nan_arg    \
        $expr_arg 
else
    echo "Directory $model_folder does not exist"
fi