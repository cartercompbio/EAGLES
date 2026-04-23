#!/bin/bash
#SBATCH --mem=32G
#SBATCH -o /cellar/users/domeyer/sbatch_outs/out/%A.%x.%a.out #STDOUT
#SBATCH -e /cellar/users/domeyer/sbatch_outs/err/%A.%x.%a.err # STDERR
#SBATCH --partition=carter-compute
#SBATCH --array=1-80%15

#parameters for scoring
scoreScript="/cellar/users/domeyer/repos/EAGLES/workflows/bin/score_other_cohort.py"
testSets="/cellar/users/domeyer/EAGLE/cohort_test_sets"
expression_base="/cellar/users/domeyer/EAGLE/cohort_expression"

modes_general=(elasticnetLDLax flipAlleleLDLax pcregressionLDLax xgbLDLax \
            elasticnetLDMed flipAlleleLDMed pcregressionLDMed xgbLDMed \
            elasticnetLDStrict flipAlleleLDStrict pcregressionLDStrict xgbLDStrict)

modes_finemapped=(elasticnetLDLax flipAlleleLDLax pcregressionLDLax xgbLDLax \
                  elasticnet flipAllele pcregression xgb)

cohorts=(gtex_eur gtex_eur_female tcga_eur tcga_eur_female)

task=$(( $SLURM_ARRAY_TASK_ID - 1 ))

if (( task < 48 )); then
    base_folder=whole_blood_nocov
    mode_idx=$(( task % 12 ))
    cohort_idx=$(( task / 12 ))
    mode=${modes_general[$mode_idx]}
else
    task2=$(( task - 48 ))
    base_folder=whole_blood_finemap
    mode_idx=$(( task2 % 8 ))
    cohort_idx=$(( task2 / 8 ))
    mode=${modes_finemapped[$mode_idx]}
fi
cohort=${cohorts[$cohort_idx]}

if [[ "$cohort" == gtex_eur* ]]; then
    expr_arg="--expr ${expression_base}/${cohort}/breast_mammary_tissue.tsv"
elif [[ "$cohort" == tcga_eur* ]]; then
    expr_arg="--expr ${expression_base}/${cohort}/breast.tsv"
else
    echo "cohort must be tcga_eur or gtex_eur. $cohort given."
    exit 1
fi


if [[ "$mode" == elasticnet* ]] || [[ "$mode" == pcregression* ]] || [[ "$mode" == flipAllele* ]]; then
    replace_nan_arg="--replace-nan 0"
else
    replace_nan_arg=""
fi

model_folder="/cellar/shared/carterlab/projects/eagle/v0.3/$base_folder/$mode/models"
score_outdir="/cellar/shared/carterlab/projects/eagle/score_tissue_from_blood/breast/${base_folder}_${mode}"
if [ -d "$model_folder" ]; then

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