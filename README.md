# EAGLES
Evaluating Aggregated Gene Level eQTL Scores

EAGLES provides pipelines to facilitate eQTL-based model training as well as scoring from pre-trained models. More detailed documentation may be found <here>. Pre-trained models are available <here>.

## Installation
```bash
git clone https://github.com/Douglas-Meyer/EAGLES.git
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
+ heritability: .tsv file with tissue,gene pairs (see EAGLES/data/<>.tsv for formatting)
+ gtexQTLfolder: folder containing single tissue eQTL tables (see EAGLES/data/<>.tsv for formatting)
    + note that file names must match tissue labels in heritability
+ gtexQTLindexFolder: folder containing .pkl index files yielded from EAGLES/workflows/bin/qtl_filter.py (index_new_eqtl_tables)
+ pfile: plink2 prefix, will be used for model training
+ train: txt file, defines the subset of samples in pfile to be used for training
+ expressionfolder: directory containing files like "gene.tsv" (see EAGLES/data/<>.tsv for formatting)
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
if [[ "$mode" == elasticnet* ]] || [[ "$mode" == pcregression* ]] || [[ "$mode" == flipAllele* ]]; then
    replace_nan_arg="--replace-nan 0"
else
    replace_nan_arg=""
fi

mkdir -p $score_outdir
conda run -n eagle python $scoreScript \
    --pgen "$testSets/$cohort.pgen" \
    --pvar "$testSets/$cohort.pvar" \
    --psam "$testSets/$cohort.psam" \
    --pindex "$testSets/$cohort.pkl" \
    --model-folder  $model_folder \
    --outdir $score_outdir     \
    $replace_nan_arg
```
