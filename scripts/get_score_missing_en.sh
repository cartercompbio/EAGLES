#!/bin/bash
#SBATCH --mem=16G
#SBATCH -o /cellar/users/domeyer/sbatch_outs/out/%A.%x.%a.out #STDOUT
#SBATCH -e /cellar/users/domeyer/sbatch_outs/err/%A.%x.%a.err # STDERR
#SBATCH --partition=carter-compute
#SBATCH --time=336:00:00
#SBATCH --array=0-21

scoreScript="/cellar/users/domeyer/repos/EAGLES/workflows/bin/score_other_cohort.py"
testSets="/cellar/users/domeyer/EAGLE/cohort_test_sets"
expression_base="/cellar/users/domeyer/EAGLE/cohort_expression"


all_basefolders=(lung lung lung thyroid thyroid thyroid thyroid thyroid thyroid breast_mammary_tissue_alphagenome breast_mammary_tissue_alphagenome breast_mammary_tissue_alphagenome breast_mammary_tissue_alphagenome breast_mammary_tissue_alphagenome breast_mammary_tissue_alphagenome breast_mammary_tissue_finemapped_all_genders breast_mammary_tissue_finemapped_all_genders breast_mammary_tissue_finemapped_all_genders colon_alphagenome colon_alphagenome colon_alphagenome colon_alphagenome)
all_modes=(elasticnetLDLax elasticnetLDLax elasticnetLDLax elasticnetLDStrict elasticnetLDStrict elasticnetLDStrict elasticnetLDLax elasticnetLDLax elasticnetLDLax elasticnet elasticnet elasticnet elasticnetLDLax elasticnetLDLax elasticnetLDLax elasticnet elasticnet elasticnet elasticnetLDLax elasticnetLDLax elasticnetLDLax elasticnetLDLax)
all_cohorts=(gtex_eur tcga_eur LUNG_CANCER gtex_eur tcga_eur THCA gtex_eur tcga_eur THCA gtex_eur tcga_eur_female BRCA gtex_eur tcga_eur_female BRCA gtex_eur tcga_eur_female BRCA gtex_eur gtex_eur tcga_eur COLORECTAL)
all_cohort_names=(gtex_eur tcga_eur LUNG_CANCER gtex_eur tcga_eur THCA gtex_eur tcga_eur THCA gtex_eur tcga_eur_female BRCA gtex_eur tcga_eur_female BRCA gtex_eur tcga_eur_female BRCA gtex_eur_transverse gtex_eur_sigmoid tcga_eur COLORECTAL)
all_exprs=(lung lung "" thyroid thyroid "" thyroid thyroid "" breast_mammary_tissue breast "" breast_mammary_tissue breast "" breast_mammary_tissue breast "" colon_transverse colon_sigmoid colorectal "")

base_folder=${all_basefolders[$SLURM_ARRAY_TASK_ID]}
mode=${all_modes[$SLURM_ARRAY_TASK_ID]}
cohort=${all_cohorts[$SLURM_ARRAY_TASK_ID]}
cohort_name=${all_cohort_names[$SLURM_ARRAY_TASK_ID]}
e=${all_exprs[$SLURM_ARRAY_TASK_ID]}


outdir="/cellar/shared/carterlab/projects/eagle/v0.3_missing_reruns/$base_folder/$mode"
model_folder="$outdir/models"   


if [[ "$cohort" == gtex* ]]; then
    expr_arg="--expr ${expression_base}/${cohort}/${e}.tsv"
    pbase="$testSets"
elif [[ "$cohort" == tcga* ]]; then
    expr_arg="--expr ${expression_base}/${cohort}/${e}.tsv"
    pbase="$testSets"
elif [[ -f "$testSets/tcga_tumor_sets/${cohort}.pkl" ]]; then
    expr_arg=""
    pbase="$testSets/tcga_tumor_sets"
else
    echo "cohort must be tcga/gtex/geuvadis or tumor. $cohort given." >&2
    exit 1
fi


if [ -d "$model_folder" ]; then
    score_outdir="$outdir/${cohort_name}"
    mkdir -p $score_outdir

    conda run -n eagle python $scoreScript \
        --pgen "$pbase/$cohort.pgen" \
        --pvar "$pbase/$cohort.pvar" \
        --psam "$pbase/$cohort.psam" \
        --pindex "$pbase/$cohort.pkl" \
        --model-folder  $model_folder \
        --outdir $score_outdir     \
        --replace-nan 0    \
        $expr_arg 
else
    echo "Directory $model_folder does not exist"
fi