// GENERAL PARAMETERS
//params.geneInfo = "/cellar/shared/carterlab/projects/eagle/grch38_gene_tss.tsv"
params.geneInfo = "/cellar/users/domeyer/EAGLE/test_expr/clean_gene_info.tsv"
params.europfile = "/carter/controlled/dbGaP/phs000178_TCGA/TOPMED_TCGA/plink2_eur_only/tcga.common.european.noimmunecancers"
params.gtexQTLfolder = "/cellar/users/domeyer/data/gtex/cis_eqtls/GTEx_EUR_slope_tables_by_ENSG_rsid_snp/by_tissue_type"
params.gtexQTLindexFolder = "/cellar/users/domeyer/data/gtex/cis_eqtls/GTEx_EUR_slope_tables_by_ENSG_rsid_snp/by_tissue_type_index"

// TCGA PARAMETERS
//params.pfile = "/carter/controlled/dbGaP/phs000178_TCGA/TOPMED_TCGA/plink2/tcga.common"
//params.expressionfolder = "/cellar/users/domeyer/data/tcga/expr_cn_by_ensg/expression"
//params.covariates  = "/cellar/users/nopopko/projects/eagles/covariates_test.csv"

// GTEX PARAMETERS
params.pfile = "/cellar/users/nopopko/projects/eagles/GTEx_plinkqc/qc_output/GTEx.qc_passed"
params.expressionfolder = "/cellar/users/domeyer/data/gtex/expression/by_tissue"
params.covariates  = "/cellar/shared/carterlab/projects/eagle/v0.2/gtex_covar/gtex.eigenvec"
params.train = "/cellar/users/domeyer/EAGLE/test_expr/eur_train_ids.txt"

params.MODE = "predixcan" //default MODE

// keep track of different runs and avoid overwriting
def timestamp = new Date().format('MMM-dd-yyyy-HH.mm')
params.outdir = "/cellar/shared/carterlab/projects/eagle/v0.2/test_out/${params.MODE}_${timestamp}"


params.predixcan = [
    window: 1000000,
    threshold: 1,
    features: "order",
    model: "elasticnet" // can change
]
params.predixcanRF = [
    window: 1000000,
    threshold: 1,
    features: "order",
    model: "rf"
]
params.magma = [
    window: 0,
    threshold: 0.999,
    features: "pca",
    model: "elasticnet" // can change
]

params.ldstrict = [
    ldWindow: "500kb",
    ldR: 0.2
]
params.ldmed = [
    ldWindow: "200kb",
    ldR: 0.5
]
params.ldlax = [
    ldWindow: "100kb",
    ldR: 0.8
]
params.ldnone = [:]

workflow {
    switch(params.MODE){
        case "predixcan":
            mode_params=params.predixcan
            mode_associated_ld = "ldnone"
            break
        case "magma":
            mode_params=params.magma
            mode_associated_ld = "ldnone"
            break
        default:
            def available_modes = ['predixcan', 'magma']
            error "Unknown MODE: '${params.MODE}'. Available modes: ${available_modes.join(', ')}"
    }
    
    // defaults to mode-associated ld parameters unless user specifies ld parameters
    def ld_mode = params.LDMODE ?: mode_associated_ld
    
    switch(ld_mode){
        case "ldlax":
            ld_params=params.ldlax
            break
        case "ldmed":
            ld_params=params.ldmed
            break
        case "ldstrict":
            ld_params=params.ldstrict
            break
        case "ldnone":
            ld_params=params.ldnone
            break
        default:
            def available_modes = ['ldlax', 'ldmed','ldstrict','ldnone']
            error "Unknown LD Mode: '${ld_mode}'. Available modes: ${available_modes.join(', ')}"
    }
    
    pfile_renamed = RENAMEVARIANTS(params.pfile)

    tissue_names = Channel
        .fromPath("${params.gtexQTLfolder}/*.tsv")
        .map { file -> file.baseName }  // Extract tissue name (remove .tsv)
        .filter { tissue ->
            // Check if corresponding .pkl file exists
            def pkl_file = new File("${params.gtexQTLindexFolder}/${tissue}.pkl")
            def expr_dir = new File("${params.expressionfolder}/${tissue}")

            // Only keep if both exist
            pkl_file.exists() && expr_dir.exists()
        }
        .take(5) // Limit to 5 tissues for now

    ginfo = Channel
        .fromPath(params.geneInfo)
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            tuple(row.ENSG, row.CHROM, row.TSS as Integer, row.STRAND, row.LENGTH as Integer)
        }
        .take(10) // Limit to 10 genes for now
        
    genes = GETSTARTSTOP(ginfo, mode_params.window) //genes: [val(ensg), val(chrom), val(start), val(end)]
    
    tissue_gene_ch = tissue_names.combine(genes) //[tis, ensg, chrom, start, end]
        .filter { tis,ensg,chrom,start,end ->
            def expr_file = new File("${params.expressionfolder}/${tis}/${ensg}.tsv")
            expr_file.exists()
        }
        
        
    switch(ld_mode){
        case "ldnone":
            variants = GETVARS(tissue_gene_ch, pfile_renamed)
            break
        default:
            variants = GETLDFILTERVARS(tissue_gene_ch, pfile_renamed, ld_params.ldWindow, ld_params.ldR)  //TODO: test this or remove it as option
    }
    
    //variants: [tis, ensg, pgen, pvar, psam]
    features = GETFEATURES(variants, mode_params.features, mode_params.threshold) // features: [tis, ensg, feats_path, loadings_path]
    
    features_for_model = features
    .map{tis,ensg,feat_path,loading_path ->
        def expression_path = file("${params.expressionfolder}/${tis}/${ensg}.tsv")
        def covariate_path = file(params.covariates)
        [tis, ensg, feat_path, expression_path, covariate_path, mode_params.model]
    }

    features_for_score = features
    .map{tis,ensg,feat_path,loading_path ->
        def model_path = file("${params.outdir}/models/${tis}_${ensg}.pkl")
        def covariate_path = file(params.covariates)
        [tis, ensg, model_path, feat_path, covariate_path]
    }

    
    models = FITMODEL(features_for_model)
    
    score_inputs = models.map { tis, ensg, model_path ->
        def covariates_path = file(params.covariates)
        def feats_path = file("${params.outdir}/features/${tis}_${ensg}_feats.tsv")
        tuple(tis, ensg, model_path, covariates_path, feats_path)
    }


    scores = MODELSCORE(score_inputs)
}

process RENAMEVARIANTS{
    cpus 1
    memory 16.GB
    
    input:
    val pfile_prefix
    
    output:
    tuple path("*.pgen"), path("*.pvar"), path("*.psam")
    
    script:
    def basename = new File(pfile_prefix).name
    """
    plink2 \\
        --pfile ${pfile_prefix} \\
        --set-all-var-ids @:#:\\\$r:\\\$a \\
        --make-pgen \\
        --new-id-max-allele-len 1000 \\
        --out ${basename}_renamed
    """
}


process GETSTARTSTOP{
    cpus 1
    memory 4.GB
    
    input:
    tuple val(ensg), val(chrom), val(tss), val(strand), val(length)
    val window
    
    output:
    tuple val(ensg), val(chrom), val(start), val(end)
    
    
    exec:
    if (params.MODE == "magma") {
        if (strand == "+") {
            start = tss - window
            end = tss + length + window
        } else if (strand == "-") {
            start = tss - length - window
            end = tss + window
        }
    } else {
        start = tss - window
        end = tss + window    
    }
    
    if (start<0){
        start = 0
    }
}


// only works if variants are named like "<chr>:<pos>:<ref>:<alt>"
// qtl_filter uses the above format to filter variants such that start <= pos <= stop
// for rsid or other naming schemes that don't use ":" this will lead to errors
process GETVARS{
    cpus 1
    memory 16.GB
    publishDir params.outdir + '/vars'
    errorStrategy 'ignore'
    
    input:
    tuple val(tis), val(ensg), val(chrom), val(start), val(stop)
    tuple path(pgen), path(pvar), path(psam)
    
    output:
    tuple val(tis), val(ensg), path("*.pgen"), path("*.pvar"), path("*.psam"), optional: true
    
    script:
    """
    python ${projectDir}/bin/qtl_filter.py \
        --qtl-folder ${params.gtexQTLfolder} \
        --index-folder ${params.gtexQTLindexFolder} \
        --gene ${ensg} \
        --tis ${tis} \
        --range ${start} ${stop} \
        --output "temp.txt"
    
    if [ -s "temp.txt" ]; then
        plink2 \
            --pfile ${pgen.baseName} \
            --chr ${chrom} \
            --from-bp ${start} \
            --to-bp ${stop} \
            --extract "temp.txt" \
            --force-intersect \
            --make-pgen \
            --out "${tis}_${ensg}"
    else
        echo "temp.txt does not exist or is empty"
    fi
    """
}

process GETLDFILTERVARS{
    cpus 1
    memory 16.GB
    publishDir params.outdir + '/filtervars'
    errorStrategy 'ignore'
    
    input:
    tuple val(ensg), val(chrom), val(start), val(stop)
    tuple path(pgen), path(pvar), path(psam)
    val(window)
    val(r2)
    
    output:
    tuple val(ensg), path("*.pgen"), path("*.pvar"), path("*.psam")
    
    script:
    """
    #!/bin/bash
    
    handle_error() {
        exit 1
    }
    
    # Set trap to catch errors
    # log genes with errors but do not emit values from this process
    trap 'handle_error' ERR
    set -e
    
    touch "temp.txt"
    python ${projectDir}/bin/qtl_filter.py \
        --qtl-folder ${params.gtexQTLfolder} \
        --index-folder ${params.gtexQTLindexFolder} \
        --gene ${ensg} \
        --tis "pantissue_no_blood_brain" \
        --output "temp.txt"
    
    plink2 \
        --pfile ${pgen.baseName} \
        --chr ${chrom} \
        --from-bp ${start} \
        --to-bp ${stop} \
        --extract "temp.txt" \
        --force-intersect \
        --make-pgen \
        --out "${ensg}_temp"

    vars=\$(mktemp)
    awk '!/^##/ && NR > 1 {print \$3}' "${ensg}_temp.pvar" > \$vars
    plink2 --indep-pairwise ${window} ${r2} --pfile ${params.europfile}  --extract \$vars --out ${ensg}
    
    plink2 --pfile "${ensg}_temp" --extract "${ensg}.prune.in" --make-pgen --out ${ensg}
    
    rm \$vars
    rm "${ensg}.prune.in"
    rm "${ensg}_temp.pvar"
    rm "${ensg}_temp.pgen"
    rm "${ensg}_temp.psam"
    rm "${ensg}_temp.log"
    
    """
}



process LDPRUNE{
    cpus 1
    memory 16.GB
    publishDir params.outdir + '/ldprune'
    
    input:
    tuple val(ensg), path(pgen), path(pvar), path(psam)
    val(window)
    val(r2)
    
    output:
    tuple val(ensg), path("*.prune.in")
    
    script:
    """
    plink2 --pfile ${params.pfile} --chr ${chrom} --from-bp ${start} --to-bp ${stop} --make-pgen --out "${ensg}_temp"
    
    vars=\$(mktemp)
    awk '!/^##/ && NR > 1 {print \$3}' "${ensg}_temp.pvar" > \$vars
    plink2 --indep-pairwise ${window} ${r2} --pfile ${params.europfile}  --extract \$vars --out ${ensg}
        
    rm \$vars
    """
}


process GETFEATURES{
    cpus 1
    memory 16.GB
    publishDir params.outdir + '/features'
    
    input:
    tuple val(tis), val(ensg), path(pgen), path(pvar), path(psam)
    val(mode)
    val(thres)
    
    output:
    tuple val(tis), val(ensg), path("${tis}_${ensg}_feats.tsv"), path("${tis}_${ensg}_loadings.tsv")
    
    script:
    """
    python ${projectDir}/bin/feature.py --pgen ${pgen} --psam ${psam} --pvar ${pvar} --method ${mode} --thres ${thres} --output "${tis}_${ensg}.tsv"
    """

}

process GETVARFEATURES{
    cpus 1
    memory 16.GB
    publishDir params.outdir + '/features'
    errorStrategy 'ignore'
    
    input:
    tuple val(tis), val(ensg), val(chrom), val(start), val(stop)
    tuple path(pgen), path(pvar), path(psam)
    val(mode)
    val(thres)
    
    output:
    tuple val(tis), val(ensg), path("${tis}_${ensg}_feats.tsv"), path("${tis}_${ensg}_loadings.tsv"), optional: true
    
    script:
    """
    
    python ${projectDir}/bin/qtl_filter.py \
        --qtl-folder ${params.gtexQTLfolder} \
        --index-folder ${params.gtexQTLindexFolder} \
        --gene ${ensg} \
        --tis ${tis} \
        --output "temp.txt"
    
    if [ -s "temp.txt" ]; then
        plink2 \
            --pfile ${pgen.baseName} \
            --chr ${chrom} \
            --from-bp ${start} \
            --to-bp ${stop} \
            --extract "temp.txt" \
            --force-intersect \
            --make-pgen \
            --out "${tis}_${ensg}"
            
        if [ -s "${tis}_${ensg}.pvar" ]; then
            python ${projectDir}/bin/feature.py \
                --pgen "${tis}_${ensg}.pgen" \
                --psam "${tis}_${ensg}.psam" \
                --pvar "${tis}_${ensg}.pvar" \
                --method ${mode} \
                --thres ${thres} \
                --output "${tis}_${ensg}"
        else
            echo "no variants"
        fi
    else
        echo "temp.txt does not exist or is empty"
    fi
    """
}

process FITMODEL{
    cpus 1
    memory 16.GB
    publishDir params.outdir + '/models'
    errorStrategy 'ignore'
    
    input:
    tuple val(tis), val(ensg), path(features), path(expression), path(covariates), val(model_type)

    output:
    tuple val(tis), val(ensg), path("*.pkl")
    
    script:
    """
    #!/bin/bash
    
    handle_error() {
        exit 1
    }
    
    # Set trap to catch errors
    # log genes with errors but do not emit values from this process
    trap 'handle_error' ERR
    set -e
    
    python ${projectDir}/bin/fit_model.py \
        --features ${features} \
        --expression ${expression} \
        --covariates ${covariates} \
        --model ${model_type} \
        --gene ${ensg} \
        --samples ${params.train} \ 
        --output "${tis}_${ensg}_${task.hash}.pkl"
    """
    
}

process MODELSCORE {
    cpus 1
    memory 32.GB
    publishDir params.outdir + '/scores'
    
    input:
    tuple val(tis), val(ensg), path(model), path(covariates), path(features)
    
    output:
    tuple val(ensg), path("*_scores.tsv")
    
    script:
    """
    python ${projectDir}/bin/model_score.py \
        --features ${features} \
        --covariates ${covariates} \
        --model ${model} \
        --output "${tis}_${ensg}_scores.tsv"
    """
}
