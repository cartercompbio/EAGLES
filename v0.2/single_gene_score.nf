params.geneInfo = "/cellar/shared/carterlab/projects/eagle/grch38_gene_tss.tsv"
params.pfile = "/carter/controlled/dbGaP/phs000178_TCGA/TOPMED_TCGA/plink2/tcga.common"
params.europfile = "/carter/controlled/dbGaP/phs000178_TCGA/TOPMED_TCGA/plink2_eur_only/tcga.common.european.noimmunecancers"
params.gtexQTLfolder = "/cellar/users/domeyer/data/gtex/cis_eqtls/GTEx_EUR_slope_tables_by_ENSG_rsid_snp/by_tissue_type"
params.gtexQTLindexFolder = "/cellar/users/domeyer/data/gtex/cis_eqtls/GTEx_EUR_slope_tables_by_ENSG_rsid_snp/by_tissue_type_index"

params.outdir="/cellar/shared/carterlab/projects/eagle/v0.2/test_out"


params.MODE = "predixcan" //default MODE

params.predixcan = [
    window: 1000000,
    threshold: 1,
    features: "order"
]

params.magma = [
    window: 0,
    threshold: 0.999,
    features: "pca"
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
        
        .take(5) // for testing purposes
        
    genes = GETSTARTSTOP(ginfo, mode_params.window)
    
    switch(ld_mode){
        case "ldnone":
            variants = GETVARS(genes)
            break
        default:
            variants = GETLDFILTERVARS(genes, ld_params.ldWindow, ld_params.ldR)
    }
                
    features = GETFEATURES(variants, mode_params.features, mode_params.threshold)
    
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
    cpus 1
    memory 16.GB
    publishDir params.outdir + '/filtervars'
    
    input:
    tuple val(ensg), val(chrom), val(start), val(stop)
    val(window)
    val(r2)
    
    output:
    tuple val(ensg), path("*.pgen"), path("*.pvar"), path("*.psam")
    
    script:
    """
    temp_file=\$(mktemp)
    python ${projectDir}/bin/qtl_filter.py \
        --qtl-folder ${gtexQTLfolder} \
        --index-folder ${gtexQTLindexFolder} \
        --gene ${ensg} \
        --tis "Colon_Sigmoid" "Colon_Transverse" \
        --output \$temp_file
    
    plink2 \
        --pfile ${params.pfile} \
        --chr ${chrom} \
        --from-bp ${start} \
        --to-bp ${stop} \
        --extract \$temp_file \
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
    
    input:
    tuple val(ensg), val(chrom), val(start), val(stop)
    
    output:
    tuple val(ensg), path("*.pgen"), path("*.pvar"), path("*.psam")
    
    script:
    """
    temp_file=\$(mktemp)
    python ${projectDir}/bin/qtl_filter.py \
        --qtl-folder ${gtexQTLfolder} \
        --index-folder ${gtexQTLindexFolder} \
        --gene ${ensg} \
        --tis "Colon_Sigmoid" "Colon_Transverse" \
        --output \$temp_file
    
    plink2 \
        --pfile ${params.pfile} \
        --chr ${chrom} \
        --from-bp ${start} \
        --to-bp ${stop} \
        --extract \$temp_file \
        --make-pgen \
        --out ${ensg}
    
    rm \$temp_file
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
    tuple val(ensg), path("*.tsv")
    
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
    tuple val(ensg), path(features), path(expression)
    path(covariates)
    val(model)
    
    output:
    tuple val(ensg), path("*.pkl")
    
    script:
    """
    python <script #1> --features ${features} --expression ${expression} --covariates ${covariates} --model <model> --output ${ensg}.pkl
    """
    //model would specify which model to use
    // e.g. rf/xgb/ridge/elasticnet...
    //
    //script will need to 
    //1. load data from features and covariates table (x)
    //2. load data from expression
    //3. fit a model of specified types
    //       -needs to support variety of common regression types
    //4. output model as pickel (.pkl) file
    
}

process MODELSCORE{
    cpus 1
    memory 16.GB
    publishDir params.outdir + '/scores'
    
    input:
    tuple val(ensg), path(features), path(model)
    path(covariates)
    
    output:
    tuple val(ensg), path("*.tsv")
    
    script:
    """
    python <script #2> --features ${features}  --covariates ${covariates} --model ${model} --output ${ensg}.tsv
    """
    //model is computed from process FITMODEL
    //
    //script will need to
    //1. load model from pkl file
    //2. compute score for each sample in features
    //       -might involve covariate values?
}
