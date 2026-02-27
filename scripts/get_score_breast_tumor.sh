#!/bin/bash
#SBATCH --mem=32G
#SBATCH -o /cellar/users/domeyer/sbatch_outs/out/%A.%x.%a.out #STDOUT
#SBATCH -e /cellar/users/domeyer/sbatch_outs/err/%A.%x.%a.err # STDERR
#SBATCH --partition=carter-compute
#SBATCH --array=1-60%15

#parameters for scoring
cohort="BRCA"

scoreScript="/cellar/users/domeyer/repos/EAGLES/workflows/bin/score_other_cohort.py"
testSets="/cellar/users/domeyer/EAGLE/cohort_test_sets/tcga_tumor_sets"

modes=(elasticnetLDLax flipAlleleLDLax pcregressionLDLax xgbLDLax randomforestLDLax elasticnetLDMed flipAlleleLDMed pcregressionLDMed xgbLDMed randomforestLDMed elasticnetLDStrict flipAlleleLDStrict pcregressionLDStrict xgbLDStrict randomforestLDStrict)
mode_idx=$(( ($SLURM_ARRAY_TASK_ID - 1) % 15 ))
mode=${modes[$mode_idx]}

base_folders=(breast_mammary_tissue breast_mammary_tissue_all_genders breast_mammary_tissue_finemapped breast_mammary_tissue_finemapped_all_genders)
base_folder_idx=$(( ($SLURM_ARRAY_TASK_ID - 1) / 15 ))
base_folder=${base_folders[$base_folder_idx]}

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
        --pgen "$testSets/$cohort.pgen" \
        --pvar "$testSets/$cohort.pvar" \
        --psam "$testSets/$cohort.psam" \
        --pindex "$testSets/$cohort.pkl" \
        --model-folder  $model_folder \
        --outdir $score_outdir     \
        $replace_nan_arg    \
else
    echo "Directory $model_folder does not exist"
fi