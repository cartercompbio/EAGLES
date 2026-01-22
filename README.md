# EAGLES
eQTLs Aggregated for Gene LEvel Scores

## Set up pixi
1. Install Pixi (if not already installed, see: https://pixi.sh/latest/installation/)
```bash
curl -fsSL https://pixi.sh/install.sh | sh
```
## run desired pipeline
2. Get pixi .lock and .toml, and nf script from github repo; .lock and .toml must be in same folder
```bash
pixi exec --manifest-path /path/to/pixi.toml nextflow run /nfpath/to/script.nf
```
