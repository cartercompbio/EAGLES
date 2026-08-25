[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21983535.svg)](https://doi.org/10.5281/zenodo.21983535)

# EAGLES
Evaluating Aggregated Gene Level eQTL Scores

EAGLES provides pipelines to facilitate eQTL-based model training as well as scoring from pre-trained models. More detailed documentation may be found <here>. Pre-trained models are available through [zenodo](https://doi.org/10.5281/zenodo.21477399).

## Installation
```bash
git clone https://github.com/cartercompbio/EAGLES.git
```

## Set up conda environment
A. Install from repository yml file (replace "eagles" with desired name)
```bash
conda env create -f EAGLES/EAGLES.yml -n eagles
```
B. create environment from scratch (replace "eagles" with desired name)
```bash
conda create -n eagles
conda activate eagles
conda install bioconda::nextflow=25.04.06 bioconda::pgenlib=0.94.1 bioconda::plink2=2.0.0a.6.9 conda-forge::python=3.13.5 conda-forge::scipy=1.16.0 conda-forge::pandas=2.3.1 conda-forge::numpy=2.3.2 conda-forge::joblib=1.5.3 conda-forge::scikit-learn=1.7.2 conda-forge::optuna=4.9.0 conda-forge::xgboost=3.3.0 conda-forge::shap=0.52.0
```

## Add EAGLES process labels to nextflow config file
see https://docs.seqera.io/nextflow/config and https://nf-co.re/docs/running/configuration/configuration-options for more details.

A. include the following in your .nextflow/config file:
```bash
includeConfig 'EAGLES/conf/EAGLES.config'
```
B. include for EAGLES-specific pipeline runs (note EAGLES will not work without process labels)
```bash
nextflow ... -c EAGLES/conf/EAGLES.config
```

## PlinkQC your genotype tables:
adapted from [Syed et al. 2025](https://github.com/meyer-lab-cshl/plinkQC)
```bash
for chr in {1..22} X; do
    nextflow run EAGLES/workflows/plinkqc.nf \
      --Cohort_name "chr${chr}" \
      --study "$pfile_to_qc" \
      --ref "$pfile_ref" \
      --qc_dir $qc_dir       
done
```
+ pfile_to_qc: plink2 genotype file (sliced from a single chromosome)
+ pfile_ref: plink2 genotype reference file (e.g. 1000 Genomes project)
+ qc_dir: output directory 


## Use EAGLES to train models
```bash
conda run -n eagles nextflow EAGLES/workflows/single_gene_scores[_train_only].nf \
    --geneInfo "$geneInfo" \
    --europfile "$europfile" \
    --heritability "$heritability" \
    --gtexQTLfolder "$gtexQTLfolder" \
    --gtexQTLindexFolder "$gtexQTLindexFolder" \
    --pfile "$pfile" \
    --train "$train" \
    --expressionfolder "$expressionfolder" \
    --mode "$mode" \
    --outdir "$outdir"
```
+ geneInfo: gene information file (see EAGLES/data/gtex_egene_gene_info.tsv for formatting)
+ europfile: plink2 prefix, will be used for LD pruning
+ heritability: .tsv file with tissue,gene pairs (see EAGLES/data/test/heritability_test.tsv for formatting)
+ gtexQTLfolder: folder containing single tissue eQTL tables (see EAGLES/data/test/gtexQTLfolder for formatting)
    + note that file names must match tissue labels in heritability (e.g. EAGLES/data/test/gtexQTLfolder/whole_blood.tsv)
+ gtexQTLindexFolder: folder containing .pkl index files yielded from EAGLES/workflows/bin/qtl_filter.py (index_new_eqtl_tables)
+ pfile: plink2 prefix, will be used for model training
+ train: txt file, defines the subset of samples in pfile to be used for training
+ expressionfolder: directory containing files like "gene.tsv" (see EAGLES/data/test/expressionfolder for formatting)
    + expected path structure is tissue/gene.tsv (e.g. EAGLES/data/test/expressionfolder/whole_blood/ENSG00000283064.tsv)
+ mode: which EAGLES mode (model type and LD strictness to apply); currently supports:
    + elasticnet
    + elasticnetLDLax
    + elasticnetLDMed
    + elasticnetLDStrict
    + flipAllele
    + flipAlleleLDLax
    + flipAlleleLDMed
    + flipAlleleLDStrict
    + pcregression
    + pcregressionLDLax
    + pcregressionLDMed
    + pcregressionLDStrict
    + xgb
    + xgbLDLax
    + xgbLDMed
    + xgbLDStrict
+ outdir: directory for pipeline outputs
    1. parameters.json: records parameters for this pipeline run
    2. vars: plink2 files sliced from input pfile used as inputs to train a model
    3. models: pkl files containing trained models
    4. scores: (not generated from single_gene_scores_train_only.nf) scores for a trained model based on full pfile cohort

## Use EAGLES to score a new cohort with previously trained models
```bash
conda run -n eagles nextflow EAGLES/workflows/score_new_cohort.nf \
    --pfile $pfile \
    --modelfolder $model_folder \
    --outdir $score_outdir \
    [ --expr $expr_table ] \
    [ --cov $covariates ] 

```
Required
+ pfile: plink2 prefix for genotypes to be scored
    + for best results, this should contain individuals from a single local ancestry label (e.g. [GRAFpop](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/GrafPop_README.html) )
+ model_folder: path to directory with trained EAGLES models
    + [24 sets of pre-trained models](https://doi.org/10.5281/zenodo.21477399)
+ score_outdir: directory where outputs will be written
    1. feature_summary.tsv: shap values representing feature importances in the scored models
    2. missing_counts.tsv: summary of model feature snps missing from input pfile
    3. model_performance.tsv: summary of model performance in input cohort
        + Y<sub>g</sub> =  W*cov + EAGLES<sub>g</sub> + ε
        + Y<sub>g</sub>: expression of gene *g*
        + EAGLES<sub>g</sub>: EAGLES prediction for gene *g*
        + cov: covariates
        + W: covariate weights
        + ε: error
    4. scores.tsv: prediction values from each model for each sample in cohort
Optional
+ covariates: tables with numerical values
    + if not provided, numeric values from pfile.psam are used
+ expr_table: tables with expression data used to evaluate model prediction accuracy
    + if not provided, model_performance.tsv will not be generated 
