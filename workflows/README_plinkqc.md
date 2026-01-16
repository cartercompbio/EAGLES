## How to run nextflow pipeline:
Run the pipeline using the following Bash loop:

```bash
for chr in {1..22} X; do
  nextflow run /carter/users/nopopko/projects/eagles/GTEx_plinkqc/plinkqc.nf \
    -c /carter/users/nopopko/projects/eagles/GTEx_plinkqc/configs/chr${chr}.config
done
```

## Final  QC data
- Passed variants: /cellar/users/nopopko/projects/eagles/GTEx_plinkqc/qc_output/GTEx_all.vpass
- Passed samples: /cellar/users/nopopko/projects/eagles/GTEx_plinkqc/qc_output/GTEx_all.spass
- pgen/pvar/psam files: /cellar/users/nopopko/projects/eagles/GTEx_plinkqc/qc_output/GTEx.qc_passed.pgen/.pvar/.psam
