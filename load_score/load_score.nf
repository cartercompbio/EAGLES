// example run: conda run -n <env> nextflow ./load_score.nf
//<env> info
//bioconda::nextflow 25.04.06
//bioconda::plink2 2.0.0a.6.9
//conda-forge::python 3.13.5
//conda-forge::scipy 1.16.0
//conda-forge::pandas 2.3.1
//conda-forge::numpy 2.3.2

//at some point reassess these and potentially shift towards input values instead
params.outdir = "/cellar/shared/carterlab/projects/eagle/test_out"
params.eqtl_tables = "/cellar/shared/carterlab/projects/eagle/eqtl_tables.txt"
params.gene_sets = "/cellar/shared/carterlab/projects/eagle/gene_sets.txt"
params.pfile = "/carter/controlled/dbGaP/phs000178_TCGA/TOPMED_TCGA/plink2/tcga.common"

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
    eqtl_geneset = eqtls
        .combine(genesets)
    
    eqtl_subset = GENESETSLICE(eqtl_geneset)
    
    ld_in = LDPRUNE(eqtl_subset)
    
    alt_al = ALTALLELE(ld_in)
    
    raw_genotype = RAWGENOTYPE(ld_in.join(alt_al))
    
    clean_genotype = CLEANGENOTYPE(raw_genotype)
    
    to_score = clean_genotype
                .join(eqtl_subset)
    
    matscores = MATMULTSCORE(to_score)
    
    directionscores = DIRECTIONSCORE(to_score)
}

//there is probably a better way to do this to reduce staging large files
//
process GENESETSLICE{
    cpus 1
    memory 16.GB
    //publishDir params.outdir + '/geneseteqtls' //probably can omit this after debugging
    
    input:
    tuple val(tissue_id), path(tissue_path, name: "tissue_input.*"), val(gene_set_id), path(gene_set_path, name: "geneset_input.*")

    output:
    tuple val("${tissue_id}.${gene_set_id}"), path("${tissue_id}.${gene_set_id}.tsv")
    
    script:
    """
    python3 ${projectDir}/bin/gene_set.py --effects ${tissue_path} --gene-set ${gene_set_path} --output "${tissue_id}.${gene_set_id}.tsv"
    """
}

process LDPRUNE{
    cpus 1
    memory 32.GB
    publishDir params.outdir + '/ld'
    
    input:
    tuple val(tissue_eqtl_id), path(tissue_eqtl_path)
    
    output:
    tuple val(tissue_eqtl_id), path("*.prune.in")
    
    script:
    """
    vars=\$(mktemp)
    awk -F '\t' '{ print \$1 }' ${tissue_eqtl_path} | tail -n +2 | sort -u > \$vars
    
    plink2 --indep-pairwise 500 0.5 --pfile ${params.pfile}  --extract \$vars --out ${tissue_eqtl_id}
    rm \$vars
    """
}

process ALTALLELE{
    cpus 1
    memory 16.GB
    publishDir params.outdir + '/ld_alt'
    
    input:
    tuple val(tissue_eqtl_id), path(tissue_eqtl_ld_in)
    
    output:
    tuple val(tissue_eqtl_id), path("${tissue_eqtl_id}.prune_alt.txt")
    
    script:
    """
    python3 ${projectDir}/bin/var_allele.py --fpath ${tissue_eqtl_ld_in} > "${tissue_eqtl_id}.prune_alt.txt"
    """
}

process RAWGENOTYPE{
    cpus 1
    memory 4.GB
    
    input:
    tuple val(tissue_eqtl_id), path(ld_in), path(alt_allele)
    
    output:
    tuple val(tissue_eqtl_id), path("${tissue_eqtl_id}.raw")
    
    script:
    """
    plink2 --export A --export-allele ${alt_allele} --pfile ${params.pfile}  --extract ${ld_in} --out ${tissue_eqtl_id}
    """
}

process CLEANGENOTYPE{
    cpus 1
    memory 4.GB
    publishDir params.outdir + '/genotype'
    
    input:
    tuple val(tissue_eqtl_id), path(raw_genotype)
    
    output:
    tuple val(tissue_eqtl_id), path("${tissue_eqtl_id}.tsv")
    
    script:
    """
    python3 ${projectDir}/bin/clean_genotype.py --inpath ${raw_genotype} --outpath "${tissue_eqtl_id}.tsv"
    """
}

process MATMULTSCORE{
    cpus 1
    memory 4.GB
    publishDir params.outdir + '/mat_mult_scores'
    
    input:
    tuple val(tissue_eqtl_id), path(clean_genotype, name: "genotype_input.*"), path(snp_gene_eqtls, name: "effect_input.*")
    
    output:
    tuple val(tissue_eqtl_id), path("${tissue_eqtl_id}.tsv")
    
    script:
    """
    python3 ${projectDir}/bin/scores.py --genotypes ${clean_genotype} --effects ${snp_gene_eqtls} --method 'matrixmult' --output "${tissue_eqtl_id}.tsv"
    
    """
}

process DIRECTIONSCORE{
    cpus 1
    memory 4.GB
    publishDir params.outdir + '/directional_scores'
    
    input:
    tuple val(tissue_eqtl_id), path(clean_genotype, name: "genotype_input.*"), path(snp_gene_eqtls, name: "effect_input.*")
    
    output:
    tuple val(tissue_eqtl_id), path("${tissue_eqtl_id}.tsv")
    
    script:
    """
    python3 ${projectDir}/bin/scores.py --genotypes ${clean_genotype} --effects ${snp_gene_eqtls} --method 'directional' --output "${tissue_eqtl_id}.tsv"
    """
}
