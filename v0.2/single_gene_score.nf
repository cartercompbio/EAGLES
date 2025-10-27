// GENERAL PARAMETERS
params.geneInfo = "/cellar/shared/carterlab/projects/eagle/grch38_gene_tss.tsv"
params.europfile = "/carter/controlled/dbGaP/phs000178_TCGA/TOPMED_TCGA/plink2_eur_only/tcga.common.european.noimmunecancers"
params.gtexQTLfolder = "/cellar/users/domeyer/EAGLE/test_expr/qtls"
params.gtexQTLindexFolder = "/cellar/users/domeyer/EAGLE/test_expr/qtls"

// TCGA PARAMETERS
//params.pfile = "/carter/controlled/dbGaP/phs000178_TCGA/TOPMED_TCGA/plink2/tcga.common"
//params.expressionfolder = "/cellar/users/domeyer/data/tcga/expr_cn_by_ensg/expression"
//params.covariates  = "/cellar/users/nopopko/projects/eagles/covariates_test.csv"

// GTEX PARAMETERS
params.pfile = "/cellar/users/nopopko/projects/eagles/GTEx_plinkqc/qc_output/GTEx.qc_passed"
params.expressionfolder = "/cellar/users/nopopko/projects/eagles/GTEx_expression"
params.covariates  = "/cellar/users/nopopko/projects/eagles/covariates_test.csv"

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

    ginfo = Channel
        .fromPath(params.geneInfo)
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            tuple(row.ENSG, row.CHROM, row.TSS as Integer, row.STRAND, row.LENGTH as Integer)
        }
        
        //.take(100) // for testing purposes
        
    genes = GETSTARTSTOP(ginfo, mode_params.window)
    
    pfile_renamed = RENAMEVARIANTS(params.pfile)
    
    switch(ld_mode){
        case "ldnone":
            variants = GETVARS(genes, pfile_renamed)
            break
        default:
            variants = GETLDFILTERVARS(genes, pfile_renamed, ld_params.ldWindow, ld_params.ldR)
    }
    
    logs = variants.map{ ensg, pgen, pvar, psam, log -> log }.collect()
    MERGELOGS(logs)
    variants_no_logs = variants.map{ ensg, pgen, pvar, psam, log -> [ensg, pgen, pvar, psam] }

                
    features = GETFEATURES(variants_no_logs, mode_params.features, mode_params.threshold)

    // feed in mode_params.model as the model type
    fitmodels = features.map { ensg, feats_path, loadings_path ->
        def expression_path = file("${params.expressionfolder}/${ensg}.tsv")
        tuple(ensg, feats_path, expression_path, mode_params.model)
    }

    models = FITMODEL(fitmodels)

    model_scores = models
        .join(features)
        .map { ensg, model_path, feats_path, loadings_path ->
            tuple(ensg, feats_path, model_path)
        }
        .set { score_inputs }

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

process GETLDFILTERVARS{
//update with the changes made to GETVARS
    cpus 1
    memory 16.GB
    publishDir params.outdir + '/filtervars'
    
    input:
    tuple val(ensg), val(chrom), val(start), val(stop)
    path(pfile_renamed)
    val(window)
    val(r2)
    
    output:
    tuple val(ensg), path("*.pgen"), path("*.pvar"), path("*.psam")
    
    script:
    """
    temp_file=\$(mktemp)
    python ${projectDir}/bin/qtl_filter.py \
        --qtl-folder ${params.gtexQTLfolder} \
        --index-folder ${params.gtexQTLindexFolder} \
        --gene ${ensg} \
        --tis "Colon_Sigmoid" "Colon_Transverse" \
        --output \$temp_file
    
    plink2 \
        --pfile ${pfile_renamed} \
        --chr ${chrom} \
        --from-bp ${start} \
        --to-bp ${stop} \
        --extract \$temp_file \
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

process GETVARS{
    cpus 1
    memory 16.GB
    publishDir params.outdir + '/vars'
    errorStrategy 'ignore'
    
    input:
    tuple val(ensg), val(chrom), val(start), val(stop)
    tuple path(pgen), path(pvar), path(psam)
    
    output:
    tuple val(ensg), path("*.pgen"), path("*.pvar"), path("*.psam"), path("*_failed.log")
    
    script:
    """
    #!/bin/bash
    
    handle_error() {
        echo "error with GETVARS for ${ensg}" > ${ensg}_failed.log
        exit 1
    }
    
    # Set trap to catch errors
    # log genes with errors but do not emit values from this process
    trap 'handle_error' ERR
    set -e
    
    # Create log file
    touch ${ensg}_failed.log
    
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
        --out ${ensg}
    
    #rm \$temp_file    
    """
}

process MERGELOGS {
    publishDir params.outdir + '/err', mode: 'copy'
    
    input:
    path logs
    
    output:
    path "combined.log"
    
    script:
    """
    cat ${logs} > combined.log
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
    tuple val(ensg), path(pgen), path(pvar), path(psam)
    val(mode)
    val(thres)
    
    output:
    tuple val(ensg), path("${ensg}_feats.tsv"), path("${ensg}_loadings.tsv")
    
    script:
    """
    python ${projectDir}/bin/feature.py --pgen ${pgen} --psam ${psam} --pvar ${pvar} --method ${mode} --thres ${thres} --output ${ensg}.tsv
    """

}

process FITMODEL{
    cpus 1
    memory 16.GB
    publishDir params.outdir + '/models'
    
    input:
    tuple val(ensg), path(features), path(expression), val(model_type)
    
    output:
    tuple val(ensg), path("*.pkl")
    
    script:
    """
    python ${projectDir}/bin/fit_model.py \
        --features ${features} \
        --expression ${expression} \
        --model ${model_type} \
        --gene ${ensg} \
        --output ${ensg}.pkl
    """
    
    //model would specify which model to use
    // e.g. rf/xgb/ridge/elasticnet...
    //
    //script will need to 
    //1. load data from features
    //2. load data from expression
    //3. fit a model of specified types
    //       -needs to support variety of common regression types
    //4. output model as pickel (.pkl) file
    
}

process MODELSCORE {
    cpus 1
    memory 32.GB
    publishDir params.outdir + '/scores'
    
    input:
    tuple val(ensg), path(features), path(model)
    
    output:
    tuple val(ensg), path("*_scores.tsv")
    
    script:
    """
    python ${projectDir}/bin/model_score.py \
        --features ${features} \
        --model ${model} \
        --output ${ensg}_scores.tsv
    """
    //model is computed from process FITMODEL
    //
    //script will need to
    //1. load model from pkl file
    //2. compute score for each sample in features
    //       -might involve covariate values?
}
