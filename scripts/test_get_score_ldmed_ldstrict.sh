#!/bin/bash
#SBATCH --mem=64G
#SBATCH -o /cellar/users/domeyer/sbatch_outs/out/%A.%x.%a.out #STDOUT
#SBATCH -e /cellar/users/domeyer/sbatch_outs/err/%A.%x.%a.err # STDERR
#SBATCH --partition=carter-compute

# Create timing log file
timing_log="/cellar/users/domeyer/EAGLE/runtime_logs/${SLURM_JOB_ID}.timing.log"

# Record start time
start_time=$(date +%s)
start_datetime=$(date '+%Y-%m-%d %H:%M:%S')
echo "Job ID: ${SLURM_JOB_ID}" > "$timing_log"
echo "Job Name: ${SLURM_JOB_NAME}" >> "$timing_log"
echo "Start Time: $start_datetime" >> "$timing_log"

#parameters for model fitting
geneInfo="/cellar/users/domeyer/EAGLE/test_expr/gtex_egene_gene_info.tsv"
europfile="/cellar/users/domeyer/EAGLE/test_expr/ld_reference/GTEx.qc_passed.EUR"
gtexQTLfolder="/cellar/users/domeyer/data/gtex/cis_eqtls/GTEx_EUR_slope_tables_by_ENSG_rsid_snp/by_tissue_type"
gtexQTLindexFolder="/cellar/users/domeyer/data/gtex/cis_eqtls/GTEx_EUR_slope_tables_by_ENSG_rsid_snp/by_tissue_type_index"
#heritability="/cellar/users/domeyer/EAGLE/test_expr/whole_blood_heritability.tsv" #debug with 
heritability="/cellar/users/domeyer/EAGLE/test_expr/whole_blood_heritability_100.tsv"
pfile="/cellar/users/nopopko/projects/eagles/GTEx_plinkqc/qc_output/GTEx.qc_passed"
expressionfolder="/cellar/users/domeyer/data/gtex/expression/by_tissue"
train="/cellar/users/domeyer/EAGLE/test_expr/eur_train_ids.txt"
#covariates="/cellar/shared/carterlab/projects/eagle/v0.2/gtex_covar/age_sex.tsv"
NFscript="/cellar/users/domeyer/repos/EAGLES/workflows/single_gene_scores_whole_blood.nf"


#parameters for scoring
scoreScript="/cellar/users/domeyer/repos/EAGLES/workflows/bin/score_other_cohort.py"
testSets="/cellar/users/domeyer/EAGLE/cohort_test_sets"
gtex_expression="/cellar/users/domeyer/EAGLE/test_expr/whole_blood_expression.tsv"
geuvadis_expression="/cellar/users/domeyer/data/1kgp/geuvadis/expression_gencode_v26/geuvadis_tpm.tsv"

for mode in elasticnetLDMed elasticnetLDStrict flipAlleleLDMed flipAlleleLDStrict pcregressionLDMed pcregressionLDStrict xgbLDMed xgbLDStrict randomforestLDMed randomforestLDStrict
do
    if [[ "$mode" == elasticnet* ]] || [[ "$mode" == pcregression* ]] || [[ "$mode" == flipAllele* ]]; then
        replace_nan_arg="--replace-nan 0"
    else
        replace_nan_arg=""
    fi
    outdir="/cellar/shared/carterlab/projects/eagle/v0.3/whole_blood_nocov/$mode"
    model_folder="$outdir/models"
    
    for cohort in geuvadis_eur geuvadis_afr gtex_eur gtex_afr
    do
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
            --model-folder  $model_folder \
            --outdir $score_outdir     \
            $replace_nan_arg    \
            --expr $expr
    done
done

# Record end time and calculate duration
end_time=$(date +%s)
end_datetime=$(date '+%Y-%m-%d %H:%M:%S')
elapsed=$((end_time - start_time))

echo "End Time: $end_datetime" >> "$timing_log"
echo "---" >> "$timing_log"
echo "Total Runtime: $(($elapsed / 3600))h $(($elapsed % 3600 / 60))m $(($elapsed % 60))s" >> "$timing_log"