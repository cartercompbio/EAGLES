#!/bin/bash
#SBATCH --mem=32G
#SBATCH -o /cellar/users/domeyer/sbatch_outs/out/%A.%x.%a.out #STDOUT
#SBATCH -e /cellar/users/domeyer/sbatch_outs/err/%A.%x.%a.err # STDERR
#SBATCH --partition=carter-compute
#SBATCH --array=1-100%20


timing_log="/cellar/users/domeyer/EAGLE/runtime_logs/get_score_breast.log"

#parameters for scoring
scoreScript="/cellar/users/domeyer/repos/EAGLES/workflows/bin/score_other_cohort.py"
testSets="/cellar/users/domeyer/EAGLE/cohort_test_sets"
expression_base="/cellar/users/domeyer/EAGLE/cohort_expression"

modes_general=(elasticnetLDLax flipAlleleLDLax pcregressionLDLax xgbLDLax \
            elasticnetLDMed flipAlleleLDMed pcregressionLDMed xgbLDMed \
            elasticnetLDStrict flipAlleleLDStrict pcregressionLDStrict xgbLDStrict)

modes_finemapped=(elasticnetLDLax flipAlleleLDLax pcregressionLDLax xgbLDLax \
                  elasticnet flipAllele pcregression xgb)

cohorts=(gtex_eur gtex_eur_female tcga_eur tcga_eur_female BRCA)

task=$(( $SLURM_ARRAY_TASK_ID - 1 ))

if (( task < 60 )); then
    base_folder1=breast_mammary_tissue_all_genders
    base_folder2=breast_mammary_tissue
    mode_idx=$(( task % 12 ))
    cohort_idx=$(( task / 12 ))
    mode=${modes_general[$mode_idx]}
else
    task2=$(( task - 60 ))
    base_folder1=breast_mammary_tissue_finemapped_all_genders
    base_folder2=breast_mammary_tissue_finemapped
    mode_idx=$(( task2 % 8 ))
    cohort_idx=$(( task2 / 8 ))
    mode=${modes_finemapped[$mode_idx]}
fi
cohort=${cohorts[$cohort_idx]}

if [[ "$cohort" == "BRCA" ]]; then
    expr_arg=""
    pbase="$testSets/tcga_tumor_sets"
elif [[ "$cohort" == gtex_eur* ]]; then
    expr_arg="--expr ${expression_base}/${cohort}/breast_mammary_tissue.tsv"
    pbase="$testSets"
elif [[ "$cohort" == tcga_eur* ]]; then
    expr_arg="--expr ${expression_base}/${cohort}/breast.tsv"
    pbase="$testSets"
else
    echo "cohort must be tcga_eur or gtex_eur or BRCA. $cohort given."
    exit 1
fi

if [[ "$mode" == elasticnet* ]] || [[ "$mode" == pcregression* ]] || [[ "$mode" == flipAllele* ]]; then
    replace_nan_arg="--replace-nan 0"
else
    replace_nan_arg=""
fi

for base_folder in $base_folder1 $base_folder2
do

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
            $expr_arg
            
        if [ -f "$score_outdir/model_performance.tsv" ]; then
            fileCreated="Success"
        else
            fileCreated="Failed"
        fi
            
        echo  -e "conda run -n eagle python $scoreScript\n    --pgen $pbase/$cohort.pgen\n    --pvar $pbase/$cohort.pvar\n    --psam $pbase/$cohort.psam\n     --pindex $pbase/$cohort.pkl\n    --model-folder  $model_folder\n    --outdir $score_outdir\n    $replace_nan_arg\n    $expr_arg\n----------------------------------------\nmodel_performance.tsv: $fileCreated\n\n" >> "$timing_log"          
                    
    else
        echo "Directory $model_folder does not exist"
    fi
    
done