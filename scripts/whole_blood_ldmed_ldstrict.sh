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

for mode in elasticnetLDMed elasticnetLDStrict flipAlleleLDMed flipAlleleLDStrict pcregressionLDMed pcregressionLDStrict xgbLDMed xgbLDStrict randomforestLDMed randomforestLDStrict
do
    outdir="/cellar/shared/carterlab/projects/eagle/v0.3/whole_blood_nocov/$mode"
    
    # Log the command
    echo "Command:" >> "$timing_log"
    echo "conda run -n eagle nextflow \"$NFscript\" \\" >> "$timing_log"
    echo "  --geneInfo \"$geneInfo\" \\" >> "$timing_log"
    echo "  --europfile \"$europfile\" \\" >> "$timing_log"
    echo "  --gtexQTLfolder \"$gtexQTLfolder\" \\" >> "$timing_log"
    echo "  --gtexQTLindexFolder \"$gtexQTLindexFolder\" \\" >> "$timing_log"
    echo "  --heritability \"$heritability\" \\" >> "$timing_log"
    echo "  --pfile \"$pfile\" \\" >> "$timing_log"
    echo "  --expressionfolder \"$expressionfolder\" \\" >> "$timing_log"
    echo "  --train \"$train\" \\" >> "$timing_log"
    echo "  --outdir \"$outdir\" \\" >> "$timing_log"
    echo "  --mode \"$mode\"" >> "$timing_log"
    echo "" >> "$timing_log"
    
    conda run -n eagle nextflow "$NFscript" \
        --geneInfo "$geneInfo" \
        --europfile "$europfile" \
        --gtexQTLfolder "$gtexQTLfolder" \
        --gtexQTLindexFolder "$gtexQTLindexFolder" \
        --heritability "$heritability" \
        --pfile "$pfile" \
        --expressionfolder "$expressionfolder" \
        --train "$train" \
        --outdir "$outdir" \
        --mode "$mode" \
        # --covariates "$covariates" \
done

# Record end time and calculate duration
end_time=$(date +%s)
end_datetime=$(date '+%Y-%m-%d %H:%M:%S')
elapsed=$((end_time - start_time))

echo "End Time: $end_datetime" >> "$timing_log"
echo "---" >> "$timing_log"
echo "Total Runtime: $(($elapsed / 3600))h $(($elapsed % 3600 / 60))m $(($elapsed % 60))s" >> "$timing_log"