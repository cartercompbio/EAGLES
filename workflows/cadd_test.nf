// GENERAL PARAMETERS
params.geneInfo = "/cellar/users/domeyer/EAGLE/test_expr/gtex_egene_gene_info.tsv"
params.europfile = "/cellar/users/domeyer/EAGLE/test_expr/ld_reference/GTEx.qc_passed.EUR"
params.gtexQTLfolder = "/cellar/users/nopopko/projects/eagles/grievous_dbs/s2f/CADD"
params.gtexQTLindexFolder = "/cellar/users/nopopko/projects/eagles/grievous_dbs/s2f/CADD"
params.heritability = "/cellar/users/domeyer/EAGLE/test_expr/whole_blood_heritability_expecto_genes.tsv"

// GTEX PARAMETERS
params.pfile = "/cellar/users/nopopko/projects/eagles/GTEx_plinkqc/qc_output/GTEx.qc_passed"
params.expressionfolder = "/cellar/users/domeyer/data/gtex/expression/by_tissue"
params.covariates = "/cellar/shared/carterlab/projects/eagle/v0.2/gtex_covar/age_sex.tsv"
params.train = "/cellar/users/domeyer/EAGLE/test_expr/eur_train_ids.txt"

params.mode = "predixcan" //default MODE

def timestamp = new Date().format('MMM-dd-yyyy-HH.mm')
params.debug = false

params.outdir = params.debug ? 
    //debug_outdirectory
    "/cellar/shared/carterlab/projects/eagle/v0.2/debug/whole_blood_cadd_test/${params.mode}_${timestamp}" :
    
    //final_outdirectory
    "/cellar/shared/carterlab/projects/eagle/v0.2/whole_blood_cadd/${params.mode}"


params.modes = [
    predixcan: [ 
        upstream: 1000000,
        downstream: 1000000,
        threshold: 1,
        model: "elasticnet",
        ldmode: "ldnone"
    ],
    
    predixcanLDStrict: [
        upstream: 1000000,
        downstream: 1000000,
        threshold: 1,
        model: "elasticnet",
        ldmode: "ldstrict"
    ],

    predixcanLDMed: [
        upstream: 1000000,
        downstream: 1000000,
        threshold: 1,
        model: "elasticnet",
        ldmode: "ldmed"
    ],
    
    predixcanLDLax: [ 
        upstream: 1000000,
        downstream: 1000000,
        threshold: 1,
        model: "elasticnet",
        ldmode: "ldlax"
    ],

    randomforestLDMed: [ 
        upstream: 1000000,
        downstream: 1000000,
        threshold: 1,
        model: "rf",
        ldmode: "ldmed"
    ],

    randomforestLDLax: [ 
        upstream: 1000000,
        downstream: 1000000,
        threshold: 1,
        model: "rf",
        ldmode: "ldlax"
    ],
    
    flipAlleleLDStrict: [
        upstream: 1000000,
        downstream: 1000000,
        threshold: 1,
        model: "flipallele",
        ldmode: "ldstrict"
    ],
    
    flipAlleleLDMed: [
        upstream: 1000000,
        downstream: 1000000,
        threshold: 1,
        model: "flipallele",
        ldmode: "ldmed"
    ],
    
    flipAlleleLDLax: [
        upstream: 1000000,
        downstream: 1000000,
        threshold: 1,
        model: "flipallele",
        ldmode: "ldlax"
    ],
    
    emagma: [
        upstream: 1000000,
        downstream: 1000000,
        threshold: 0.999,
        model: "pcr",
        ldmode: "ldnone"
    ],
    
    emagmaLDStrict: [
        upstream: 1000000,
        downstream: 1000000,
        threshold: 0.999,
        model: "pcr",
        ldmode: "ldstrict"
    ],
    
    emagmaLDMed: [
        upstream: 1000000,
        downstream: 1000000,
        threshold: 0.999,
        model: "pcr",
        ldmode: "ldmed"
    ],
    
    emagmaLDLax: [
        upstream: 1000000,
        downstream: 1000000,
        threshold: 0.999,
        model: "pcr",
        ldmode: "ldlax"
    ],
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

include { RENAMEVARIANTS; GETSTARTSTOP; GETVARS_BATCH } from '../modules/preprocessing'
include { FITMODEL; MODELSCORE } from '../modules/modeling'


workflow {
    if (!params.modes.containsKey(params.mode)) {
        error "Unknown mode: '${params.mode}'. Available modes: ${params.modes.keySet().join(', ')}"
    }

    mode_params = params.modes[params.mode]
    
    switch(mode_params.ldmode){
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
            error "Unknown LD Mode: '${mode_params.ldmode}'. Available modes: ${available_modes.join(', ')}"
    }
    
    pfile_renamed = RENAMEVARIANTS(params.pfile)

    tissue_names = Channel.of("whole_blood")

    heritable_gene = Channel
        .fromPath(params.heritability)
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            row.ENSG
        }
        
    gene_info_ch = Channel.fromPath(params.geneInfo)
    
    GETSTARTSTOP(
        gene_info_ch,
        mode_params.upstream,
        mode_params.downstream
    )
    
    genes = GETSTARTSTOP.out
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            tuple(row.ENSG, row.CHROM, row.start, row.end)
        }
        .join(heritable_gene)        
    
    tissue_gene_ch = tissue_names.combine(genes) //[tis, ensg, chrom, start, end]
        .filter { tis,ensg,chrom,start,end ->
            def expr_file = new File("${params.expressionfolder}/${tis}/${ensg}.tsv")
            expr_file.exists()
        }

    heritable_tis_gene = Channel
        .fromPath(params.heritability)
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            tuple(row.TISSUE, row.ENSG)
        }
    
    if (params.debug) {
        filtered_tis_gene_ch = tissue_gene_ch.join(heritable_tis_gene, by: [0,1])//[tis, ensg, chrom, start, end]
            .groupTuple(by: [1, 2, 3, 4]) 
            .take(30)  
    } else {
        filtered_tis_gene_ch = tissue_gene_ch.join(heritable_tis_gene, by: [0,1])//[tis, ensg, chrom, start, end]
            .groupTuple(by: [1, 2, 3, 4]) 
    }
    
        
    switch(mode_params.ldmode){
        case "ldnone":
            variant_res = GETVARS_BATCH(filtered_tis_gene_ch, pfile_renamed, false, "", "")
            break
        default:
            variant_res = GETVARS_BATCH(filtered_tis_gene_ch, pfile_renamed, true,  ld_params.ldWindow, ld_params.ldR)
            break
    }
    
    variants = variant_res
        .transpose()
        .map { pgen, psam, pvar ->
            def filename = pgen.baseName  
            def match = filename =~ /^(.+)_(ENSG\d+)$/
            def tis = match[0][1] 
            def ensg = match[0][2]  
            [tis, ensg, pgen, psam, pvar]
        }
        

    variants_for_model = variants
        .map { tis, ensg, pgen, psam, pvar ->
            def expression_path = file("${params.expressionfolder}/${tis}/${ensg}.tsv", checkIfExists: true)
            def covariate_path = file(params.covariates, checkIfExists: true)
            def qtl_path = file("${params.gtexQTLfolder}/whole_blood.tsv", checkIfExists: true)
            def qtl_index_path = file("${params.gtexQTLindexFolder}/whole_blood.pkl", checkIfExists: true)
    
            [tis, ensg, file(pgen, checkIfExists: true), file(psam, checkIfExists: true), file(pvar, checkIfExists: true),
             expression_path, covariate_path, qtl_path, qtl_index_path]
        }
        .filter { tis, ensg, pgen, psam, pvar, expr, covar, qtl, qtl_index ->
            def all_exist = pgen.exists() && psam.exists() && pvar.exists() &&
                            expr.exists() && covar.exists() && qtl.exists() && qtl_index.exists()
            if (!all_exist) {
                println "Skipping $ensg in $tis — missing some input files"
            }
            all_exist
    }




    models = FITMODEL(variants_for_model, mode_params.model, mode_params.threshold, params.train)
    
    variants_for_scores = variants
        .join(models, by: [0,1])
        .map{tis, ensg, pgen, psam, pvar, model_path ->
        def covariate_path = file(params.covariates)
        [tis, ensg, pgen,psam,pvar, model_path, covariate_path]
    } 
    scores = MODELSCORE(variants_for_scores)
}
