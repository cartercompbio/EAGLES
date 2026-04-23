#!/bin/bash
#SBATCH --mem=32G
#SBATCH -o /cellar/users/domeyer/sbatch_outs/out/%A.%x.%a.out #STDOUT
#SBATCH -e /cellar/users/domeyer/sbatch_outs/err/%A.%x.%a.err # STDERR
#SBATCH --partition=carter-compute
#SBATCH --array=1-80%15

#parameters for scoring
scoreScript="/cellar/users/domeyer/repos/EAGLES/workflows/bin/score_other_cohort.py"
testSets="/cellar/users/domeyer/EAGLE/cohort_test_sets"
gtex_expression="/cellar/users/domeyer/EAGLE/test_expr/whole_blood_expression.tsv"
geuvadis_expression="/cellar/users/domeyer/data/1kgp/geuvadis/expression_gencode_v26/geuvadis_tpm.tsv"

cohorts=(geuvadis_eur geuvadis_afr gtex_eur gtex_afr)
modes_general=(elasticnetLDLax flipAlleleLDLax pcregressionLDLax xgbLDLax \
            elasticnetLDMed flipAlleleLDMed pcregressionLDMed xgbLDMed \
            elasticnetLDStrict flipAlleleLDStrict pcregressionLDStrict xgbLDStrict)

modes_finemapped=(elasticnetLDLax flipAlleleLDLax pcregressionLDLax xgbLDLax \
                  elasticnet flipAllele pcregression xgb)

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

if [[ "$mode" == elasticnet* ]] || [[ "$mode" == pcregression* ]] || [[ "$mode" == flipAllele* ]]; then
    replace_nan_arg="--replace-nan 0"
else
    replace_nan_arg=""
fi

if [[ "$cohort" == geuvadis* ]]; then
    expr=$geuvadis_expression
elif [[ "$cohort" == gtex* ]]; then
    expr=$gtex_expression
else
    echo "Error: cohort must start with 'geuvadis' or 'gtex', got: $cohort" >&2
    exit 1
fi

outdir="/cellar/shared/carterlab/projects/eagle/v0.3/$base_folder/$mode"
model_folder="$outdir/models"
score_outdir="$outdir/${cohort}"

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
        --expr $expr
else
    echo "Directory $model_folder does not exist"
fi

