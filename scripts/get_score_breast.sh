#!/bin/bash
#SBATCH --mem=32G
#SBATCH -o /cellar/users/domeyer/sbatch_outs/out/%A.%x.%a.out #STDOUT
#SBATCH -e /cellar/users/domeyer/sbatch_outs/err/%A.%x.%a.err # STDERR
#SBATCH --partition=carter-compute
#SBATCH --array=1-45%10

#parameters for scoring
scoreScript="/cellar/users/domeyer/repos/EAGLES/workflows/bin/score_other_cohort.py"
testSets="/cellar/users/domeyer/EAGLE/cohort_test_sets"
gtex_expression="/cellar/users/domeyer/EAGLE/cohort_expression/gtex_eur/breast_mammary_tissue.tsv"
tcga_expression="/cellar/users/domeyer/data/tcga/EAGLES/expression_by_tissue_type/breast.tsv"

cohorts=(tcga_eur_female gtex_eur_female gtex_eur)
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

if [[ "$cohort" == tcga* ]]; then
    expr=$tcga_expression
elif [[ "$cohort" == gtex* ]]; then
    expr=$gtex_expression
else
    echo "Error: cohort must start with 'tcga' or 'gtex', got: $cohort" >&2
    exit 1
fi


for base_folder in breast_mammary_tissue breast_mammary_tissue_all_genders
do

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
    
done