#!/bin/bash
#SBATCH --mem=16G
#SBATCH -o /cellar/users/domeyer/sbatch_outs/out/%A.%x.%a.out #STDOUT
#SBATCH -e /cellar/users/domeyer/sbatch_outs/err/%A.%x.%a.err # STDERR
#SBATCH --partition=carter-compute
#SBATCH --time=336:00:00
#SBATCH --array=1-16

scoreScript="/cellar/users/domeyer/repos/EAGLES/workflows/bin/score_other_cohort.py"
testSets="/cellar/users/domeyer/EAGLE/cohort_test_sets"
expression_base="/cellar/users/domeyer/EAGLE/cohort_expression"

modes=(elasticnetLDStrict flipAlleleLDStrict pcregressionLDStrict xgbLDStrict)
basefolders=(skin_sun_exposed_lower_leg_218 whole_blood_218 skin_sun_exposed_lower_leg_889 whole_blood_889)
task=$(( $SLURM_ARRAY_TASK_ID - 1 ))
mode_idx=$(( task % 4 ))
mode=${modes[$mode_idx]}
basefolder_idx=$(( task / 4 ))
basefolder=${basefolders[$basefolder_idx]}

cohort="icb"

if [[ "$basefolder" == whole_blood* ]]; then
    expr_tis="whole_blood"
elif [[ "$basefolder" == skin* ]]; then
    expr_tis="skin_sun_exposed_lower_leg"
else
    echo "tis must be whole_blood, or skin. $basefolder given." >&2
    exit 1
fi

if [[ "$cohort" == "icb" ]]; then
    expr_arg=""
else
    expr_arg="--expr ${expression_base}/${cohort}/${expr_tis}.tsv"
fi
    
outdir="/cellar/shared/carterlab/projects/eagle/icb_scores/$basefolder/$mode"
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
        --pgen "$testSets/$cohort.pgen" \
        --pvar "$testSets/$cohort.pvar" \
        --psam "$testSets/$cohort.psam" \
        --pindex "$testSets/$cohort.pkl" \
        --model-folder  $model_folder \
        --outdir $score_outdir     \
        $replace_nan_arg    \
        $expr_arg 
else
    echo "Directory $model_folder does not exist"
fi

