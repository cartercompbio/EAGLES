workflow {
    eqtls = Channel
        .fromPath(params.eqtl_tables)
        .splitText()
        .map { line ->
            f = file(line.strip())
            [f.name[0..-5], f]  // remove .tsv to get tissue_id
        }
        
    genesets = Channel
        .fromPath(params.gene_sets)
        .splitText()
        .map {line ->
            f = file(line.strip())
            [f.name[0..-5], f] //remove .txt to get gene_set_id
        }
        
    (ld_in, ld_out) = LDPRUNE(eqtls)
    (ld_bed, ld_bim, ld_fam) = BED(eqtls.join(ld_in))
    
     
}

process LDPRUNE{
    cpus 1
    memory 32.GB
    publishDir params.outdir + '/ld'
    
    input:
    tuple val(tissue_id), path(tissue_path)
    
    output:
    tuple val(tissue_id), path("*.prune.in")
    path "${tissue_id}.prune.out"
    
    script:
    """
    vars=\$(mktemp)
    awk -F '\t' '{ print \$2 }' ${tissue_path} | tail -n +2 | sort -u > \$vars
    
    plink2 --indep-pairwise 500 0.5 --pfile ${params.pfile}  --extract \$vars --out ${tissue_id}
    rm \$vars
    """
}

process BED{
    cpus 1
    memory 32.GB
    publishDir '/cellar/users/domeyer/restricted/temp'
    
    input:
    tuple val(tissue_id), path(tissue_path), path(ld_in)
    
    output:
    path "${tissue_id}.bed"
    path "${tissue_id}.bim"
    path "${tissue_id}.fam"
    
    script:
    """
    plink2 --make-bed --pfile ${params.pfile}  --extract ${ld_in} --out ${tissue_id}
    
    """
}
