#!/bin/bash
#SBATCH --mem=32G
#SBATCH -o /cellar/users/domeyer/sbatch_outs/out/%A.%x.%a.out #STDOUT
#SBATCH -e /cellar/users/domeyer/sbatch_outs/err/%A.%x.%a.err # STDERR
#SBATCH --partition=carter-compute
#SBATCH --array=1-75%15

#parameters for scoring

scoreScript="/cellar/users/domeyer/repos/EAGLES/workflows/bin/score_other_cohort.py"
testSets="/cellar/users/domeyer/EAGLE/cohort_test_sets"

gtex_expression="/cellar/users/domeyer/EAGLE/cohort_expression/gtex_eur/lung.tsv"
tcga_expression="/cellar/users/domeyer/data/tcga/EAGLES/expression_by_tissue_type/lung.tsv"

modes_general=(elasticnetLDLax flipAlleleLDLax pcregressionLDLax xgbLDLax randomforestLDLax \
            elasticnetLDMed flipAlleleLDMed pcregressionLDMed xgbLDMed randomforestLDMed \
            elasticnetLDStrict flipAlleleLDStrict pcregressionLDStrict xgbLDStrict randomforestLDStrict)

modes_finemapped=(elasticnetLDLax flipAlleleLDLax pcregressionLDLax xgbLDLax randomforestLDLax \
                  elasticnet flipAllele pcregression xgb randomforest)

cohorts=(gtex_eur tcga_eur LUNG_CANCER)

task=$(( $SLURM_ARRAY_TASK_ID - 1 ))

if (( task < 45 )); then
    base_folder=lung
    mode_idx=$(( task % 15 ))
    cohort_idx=$(( task / 15 ))
    mode=${modes_general[$mode_idx]}
else
    task2=$(( task - 45 ))
    base_folder=lung_finemapped
    mode_idx=$(( task2 % 10 ))
    cohort_idx=$(( task2 / 10 ))
    mode=${modes_finemapped[$mode_idx]}
fi
cohort=${cohorts[$cohort_idx]}



if [[ "$cohort" == LUNG_CANCER]]; then
    expr_arg=""
    pbase="$testSets/tcga_tumor_sets"
elif [["$cohort" == gtex*]]; then
    expr_arg="--expr $gtex_expression"
    pbase="$testSets"
elif [["$cohort" == tcga*]]; then
    expr_arg="--expr $tcga_expression"
    pbase="$testSets"
else
    echo "cohort must be tcga or gtex or LUNG_CANCER. $cohort given."
    exit 1
fi


if [[ "$mode" == elasticnet* ]] || [[ "$mode" == pcregression* ]] || [[ "$mode" == flipAllele* ]]; then
    replace_nan_arg="--replace-nan 0"
else
    replace_nan_arg=""
fi


outdir="/cellar/shared/carterlab/projects/eagle/v0.3/$base_folder/$mode"
model_folder="$outdir/models"
score_outdir="$outdir/${cohort}"
if [ -d "$model_folder" ]; then

    mkdir -p $score_outdir

    conda run -n eagle python $scoreScript \
        --pgen "$pbase/$cohort.pgen" \
        --pvar "$pbase/$cohort.pvar" \
        --psam "$pbase/$cohort.psam" \
        --pindex "$pbase/$cohort.pkl" \
        --model-folder  $model_folder \
        --outdir $score_outdir     \
        $replace_nan_arg    \
        $expr_arg\
        
else
    echo "Directory $model_folder does not exist"
fi