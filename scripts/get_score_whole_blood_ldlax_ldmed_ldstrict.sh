#!/bin/bash
#SBATCH --mem=32G
#SBATCH -o /cellar/users/domeyer/sbatch_outs/out/%A.%x.%a.out #STDOUT
#SBATCH -e /cellar/users/domeyer/sbatch_outs/err/%A.%x.%a.err # STDERR
#SBATCH --partition=carter-compute
#SBATCH --array=1-60%5

#parameters for scoring
scoreScript="/cellar/users/domeyer/repos/EAGLES/workflows/bin/score_other_cohort.py"
testSets="/cellar/users/domeyer/EAGLE/cohort_test_sets"
gtex_expression="/cellar/users/domeyer/EAGLE/test_expr/whole_blood_expression.tsv"
geuvadis_expression="/cellar/users/domeyer/data/1kgp/geuvadis/expression_gencode_v26/geuvadis_tpm.tsv"

cohorts=(geuvadis_eur geuvadis_afr gtex_eur gtex_afr)
modes=(elasticnetLDLax flipAlleleLDLax pcregressionLDLax xgbLDLax randomforestLDLax elasticnetLDMed flipAlleleLDMed pcregressionLDMed xgbLDMed randomforestLDMed elasticnetLDStrict flipAlleleLDStrict pcregressionLDStrict xgbLDStrict randomforestLDStrict)

cohort_idx=$(( ($SLURM_ARRAY_TASK_ID - 1) / 15 ))
mode_idx=$(( ($SLURM_ARRAY_TASK_ID - 1) % 15 ))

cohort=${cohorts[$cohort_idx]}
mode=${modes[$mode_idx]}

if [[ "$mode" == elasticnet* ]] || [[ "$mode" == pcregression* ]] || [[ "$mode" == flipAllele* ]]; then
    replace_nan_arg="--replace-nan 0"
else
    replace_nan_arg=""
fi

outdir="/cellar/shared/carterlab/projects/eagle/v0.3/whole_blood_nocov/$mode"
model_folder="$outdir/models"

score_outdir="$outdir/${cohort}"
mkdir -p $score_outdir


if [[ "$cohort" == geuvadis* ]]; then
    expr=$geuvadis_expression
elif [[ "$cohort" == gtex* ]]; then
    expr=$gtex_expression
else
    echo "Error: cohort must start with 'geuvadis' or 'gtex', got: $cohort" >&2
    exit 1
fi

conda run -n eagle python $scoreScript \
    --pgen "$testSets/$cohort.pgen" \
    --pvar "$testSets/$cohort.pvar" \
    --psam "$testSets/$cohort.psam" \
    --pindex "$testSets/$cohort.pkl" \
    --model-folder  $model_folder \
    --outdir $score_outdir     \
    $replace_nan_arg    \
    --expr $expr