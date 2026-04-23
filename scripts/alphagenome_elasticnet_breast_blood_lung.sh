#!/bin/bash
#SBATCH --mem=64G
#SBATCH -o /cellar/users/domeyer/sbatch_outs/out/%A.%x.%a.out #STDOUT
#SBATCH -e /cellar/users/domeyer/sbatch_outs/err/%A.%x.%a.err # STDERR
#SBATCH --partition=carter-compute
#SBATCH --time=336:00:00

# Create timing log file
timing_log="/cellar/users/domeyer/EAGLE/runtime_logs/${SLURM_JOB_ID}.timing.log"

# Record start time
start_time=$(date +%s)
start_datetime=$(date '+%Y-%m-%d %H:%M:%S')
echo "Job ID: ${SLURM_JOB_ID}" > "$timing_log"
echo "Job Name: ${SLURM_JOB_NAME}" >> "$timing_log"
echo "Start Time: $start_datetime" >> "$timing_log"

#parameters for pipeline
geneInfo="/cellar/users/domeyer/EAGLE/test_expr/gtex_egene_gene_info.tsv"
europfile="/cellar/users/domeyer/EAGLE/test_expr/ld_reference/GTEx.qc_passed.EUR"
gtexQTLfolder="/cellar/users/domeyer/data/gtex/alphagenome_high_quantile"
gtexQTLindexFolder="/cellar/users/domeyer/data/gtex/alphagenome_high_quantile"

pfile="/cellar/users/nopopko/projects/eagles/GTEx_plinkqc/qc_output/GTEx.qc_passed"
expressionfolder="/cellar/users/domeyer/data/gtex/expression/by_tissue"
train="/cellar/users/domeyer/EAGLE/test_expr/eur_train_ids.txt"
NFscript="/cellar/users/domeyer/repos/EAGLES/workflows/single_gene_scores.nf"

for tissue in whole_blood breast_mammary_tissue lung
do
    

    if [[ "$tissue" == "whole_blood" ]]; then
        heritability="/cellar/users/domeyer/EAGLE/test_expr/whole_blood_heritability.tsv"
    elif [[ "$tissue" == "breast_mammary_tissue" ]]; then
        heritability="/cellar/users/domeyer/EAGLE/heritability_tables_solid_tissues/breast_mammary_tissue_heritability.tsv"
    elif [[ "$tissue" == "lung" ]]; then
        heritability="/cellar/users/domeyer/EAGLE/heritability_tables_solid_tissues/lung_heritability.tsv"
    else
        echo "tissue must be whole_blood, lung or breast_mammary_tissue. $tissue given."
        exit 1
    fi
    
    for mode in elasticnet elasticnetLDLax elasticnetLDMed elasticnetLDStrict
    do
    
        outdir="/cellar/shared/carterlab/projects/eagle/v0.3/${tissue}_alphagenome/$mode"

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
            --mode "$mode"
    done
done

# Record end time and calculate duration
end_time=$(date +%s)
end_datetime=$(date '+%Y-%m-%d %H:%M:%S')
elapsed=$((end_time - start_time))

echo "End Time: $end_datetime" >> "$timing_log"
echo "---" >> "$timing_log"
echo "Total Runtime: $(($elapsed / 3600))h $(($elapsed % 3600 / 60))m $(($elapsed % 60))s" >> "$timing_log"