#!/bin/bash

#SBATCH --mem=32G
#SBATCH -o /cellar/users/domeyer/sbatch_outs/out/%A.%x.%a.out #STDOUT
#SBATCH -e /cellar/users/domeyer/sbatch_outs/err/%A.%x.%a.err # STDERR
#SBATCH --partition=carter-compute
#SBATCH --array=1-2

#cohorts=(gtex_train gtex_all_eur gtex_eur gtex_afr geuvadis_eur geuvadis_afr)
#cohort_idx=$(( ($SLURM_ARRAY_TASK_ID - 1) % 6 ))
#cohort=${cohorts[$cohort_idx]}
cohort=gtex_all_eur

modes=(eqtl finemap)
#mode_idx=$(( ($SLURM_ARRAY_TASK_ID - 1) / 6 ))
mode_idx=$(( ($SLURM_ARRAY_TASK_ID - 1)))
mode=${modes[$mode_idx]}


gtex_expression="/cellar/users/domeyer/EAGLE/test_expr/whole_blood_expression.tsv"
geuvadis_expression="/cellar/users/domeyer/data/1kgp/geuvadis/expression_gencode_v26/geuvadis_tpm.tsv"

if [[ "$cohort" == geuvadis* ]]; then
    expr=$geuvadis_expression
    testSets="/cellar/users/domeyer/EAGLE/cohort_test_sets"
elif [[ "$cohort" == gtex_train ]]; then
    expr=$gtex_expression
    testSets="/cellar/users/domeyer/EAGLE/cohort_test_sets/gtex_train"
elif [[ "$cohort" == gtex_all_eur ]]; then
    expr=$gtex_expression
    testSets="/cellar/users/domeyer/EAGLE/cohort_test_sets/gtex_all_eur"
elif [[ "$cohort" == gtex* ]]; then 
    expr=$gtex_expression
    testSets="/cellar/users/domeyer/EAGLE/cohort_test_sets"
else
    echo "Error: cohort must start with 'geuvadis' or 'gtex', got: $cohort" >&2
    exit 1
fi

if [[ "$mode" == finemap ]]; then
    eqtlFolder="/cellar/users/domeyer/data/gtex/finemapped_v10/by_tissue"
else
    eqtlFolder="/cellar/users/domeyer/data/gtex/cis_eqtls/GTEx_EUR_slope_tables_by_ENSG_rsid_snp/by_tissue_type"
fi

scoreScript="/cellar/users/domeyer/repos/EAGLES/workflows/bin/compare_qtl_cohorts.py"

conda run -n eagle python $scoreScript \
    --pgen "$testSets/$cohort.pgen" \
    --pvar "$testSets/$cohort.pvar" \
    --psam "$testSets/$cohort.psam" \
    --pindex "$testSets/$cohort.pkl" \
    --expr $expr \
    --eqtl-folder $eqtlFolder \
    --output "/cellar/users/domeyer/EAGLE/Figures/rough_tables/${mode}_${cohort}.tsv"