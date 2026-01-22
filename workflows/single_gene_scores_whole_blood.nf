// GENERAL PARAMETERS
params.geneInfo = "/cellar/users/domeyer/EAGLE/test_expr/gtex_egene_gene_info.tsv"
params.europfile = "/cellar/users/domeyer/EAGLE/test_expr/ld_reference/GTEx.qc_passed.EUR"
params.gtexQTLfolder = "/cellar/users/domeyer/data/gtex/cis_eqtls/GTEx_EUR_slope_tables_by_ENSG_rsid_snp/by_tissue_type"
params.gtexQTLindexFolder = "/cellar/users/domeyer/data/gtex/cis_eqtls/GTEx_EUR_slope_tables_by_ENSG_rsid_snp/by_tissue_type_index"
//


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
    "/cellar/shared/carterlab/projects/eagle/v0.3/debug/${params.mode}_${timestamp}" :
    
    //final_outdirectory
    "/cellar/shared/carterlab/projects/eagle/v0.3/whole_blood_test100/${params.mode}"
    
    
params.heritability = params.debug ? 
    //debug file
    "/cellar/users/domeyer/EAGLE/test_expr/whole_blood_heritability_100.tsv" :
    
    //full file
    "/cellar/users/domeyer/EAGLE/test_expr/whole_blood_heritability.tsv"

params.modes = [
    elasticnet: [ 
        upstream: 1000000,
        downstream: 1000000,
        threshold: 1,
        model: "elasticnet",
        ldmode: "ldnone"
    ],
    
    elasticnetLDStrict: [
        upstream: 1000000,
        downstream: 1000000,
        threshold: 1,
        model: "elasticnet",
        ldmode: "ldstrict"
    ],

    elasticnetLDMed: [
        upstream: 1000000,
        downstream: 1000000,
        threshold: 1,
        model: "elasticnet",
        ldmode: "ldmed"
    ],
    
    elasticnetLDLax: [ 
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
    
    xgbLDMed: [ 
        upstream: 1000000,
        downstream: 1000000,
        threshold: 1,
        model: "xgb",
        ldmode: "ldmed"
    ],

    xgbLDLax: [ 
        upstream: 1000000,
        downstream: 1000000,
        threshold: 1,
        model: "xgb",
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
    
    pcregression: [
        upstream: 1000000,
        downstream: 1000000,
        threshold: 0.999,
        model: "pcr",
        ldmode: "ldnone"
    ],
    
    pcregressionLDStrict: [
        upstream: 1000000,
        downstream: 1000000,
        threshold: 0.999,
        model: "pcr",
        ldmode: "ldstrict"
    ],
    
    pcregressionLDMed: [
        upstream: 1000000,
        downstream: 1000000,
        threshold: 0.999,
        model: "pcr",
        ldmode: "ldmed"
    ],
    
    pcregressionLDLax: [
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
    

    filtered_tis_gene_ch = tissue_gene_ch.join(heritable_tis_gene, by: [0,1])//[tis, ensg, chrom, start, end]
        .groupTuple(by: [1, 2, 3, 4]) 
    
        
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
        .map{tis,ensg,pgen,psam,pvar ->
            def expression_path = file("${params.expressionfolder}/${tis}/${ensg}.tsv")
            def covariate_path = file(params.covariates)
            def qtl_path = file("${params.gtexQTLfolder}/${tis}.tsv")
            [tis,ensg,pgen,psam,pvar,expression_path,covariate_path, qtl_path]
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
