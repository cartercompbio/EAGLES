# EAGLES
Evaluating Aggregated Gene Level eQTL Scores

## Set up pixi
1. Install Pixi (if not already installed, see: https://pixi.sh/latest/installation/)
```bash
curl -fsSL https://pixi.sh/install.sh | sh
```
2. Create environment from root of respository (where pixi.toml/.lock live)
```bash
pixi install
```
3. Run the nextflow pipeline
```bash
pixi run score
```

## to replace 2+3 probably
2. Get pixi .lock and .toml, and nf script from github repo; .lock and .toml must be in same folder
```bash
pixi exec --manifest-path /path/to/pixi.toml nextflow run /nfpath/to/script.nf
```

## Set up conda environment with necessary packages
conda install bioconda::nextflow=25.04.06 bioconda::plink2=2.0.0a.6.9 conda-forge::python=3.13.5 conda-forge::scipy=1.16.0 conda-forge::pandas=2.3.1 conda-forge::numpy=2.3.2 conda-forge::joblib=1.5.3 conda-forge::scikit-learn=1.9.0 conda-forge::optuna=4.9.0 conda-forge::xgboost=3.3.0 conda-forge::shap=0.52.0

## Set up nextflow config file
1. touch /cellar/users/<**yours**>/.nextflow/config
2. vim /cellar/users/<**yours**>/.nextflow/config and then paste in the following (changing to yours)
3. esc then ":wq" to save config file

resume = true
dag.overwrite = true
workDir='/cellar/users/domeyer/nextflow/work' **change to yours**

process {
    cpus = 1
    memory = '6GB'
    executor = 'slurm'
    queue = 'carter-compute'
    errorStrategy = 'finish'
}

executor {
    name = 'slurm'
    queueSize = 32
}

notification {
    enabled = true 
    to = 'domeyer@ucsd.edu' **change to yours**
}

## modify wrapper bash script (run_load_score.sh)
#! /bin/bash\
#SBATCH --mem=64G\
#SBATCH -o /cellar/users/domeyer/sbatch_outs/out/%A.%x.%a.out #STDOUT\
#SBATCH -e /cellar/users/domeyer/sbatch_outs/err/%A.%x.%a.err # STDERR\
#SBATCH --partition=carter-compute\
conda run -n eagle nextflow /cellar/users/domeyer/EAGLE/load_score/load_score.nf\
\
should become
\
#! /bin/bash\
#SBATCH --mem=64G\
#SBATCH -o <**output_folder**>\
#SBATCH -e <**error_folder**>\
#SBATCH --partition=carter-compute\
conda run -n <**conda_env**> nextflow <**/path/to/load_score.nf**>

## change output directory (to not overwrite existing results)
in load_score.nf change line 11 change
+ params.outdir = "/cellar/shared/carterlab/projects/eagle/test_out"\
to
+ params.outdir = "/cellar/shared/carterlab/projects/eagle/test_out_noa"

## run pipeline
should be as simple as sbatch /path/to/run_load_score.sh
