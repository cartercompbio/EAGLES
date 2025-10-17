## How to run nextflow pipeline:
Run the pipeline using the following Bash loop:

```bash
for chr in {1..22}; do
  nextflow run /carter/users/nopopko/projects/eagles/GTEx_plinkqc/plinkqc.nf \
    -c /carter/users/nopopko/projects/eagles/GTEx_plinkqc/configs/chr${chr}.config
done
```

## Final  QC data
- Passed variants: GTEx_all.vpass
- Passed samples: GTEx_all.spass
- pgen/pvar/psam files: GTEx.qc_passed.pgen/.pvar/.psam
