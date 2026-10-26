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
 
    xgb: [ 
        upstream: 1000000,
        downstream: 1000000,
        threshold: 1,
        model: "xgb",
        ldmode: "ldnone"
    ],
    xgbLDStrict: [ 
        upstream: 1000000,
        downstream: 1000000,
        threshold: 1,
        model: "xgb",
        ldmode: "ldstrict"
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
    
    flipAllele: [ 
        upstream: 1000000,
        downstream: 1000000,
        threshold: 1,
        model: "flipallele",
        ldmode: "ldnone"
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


include { RENAMEVARIANTS; MAFFILTER; PREPROCESSEXPR; GENELISTPICKLE;  GETSTARTSTOP; GETVARS; GETQTLPICKLE; } from '../modules/preprocessing'
include { FITMODEL } from '../modules/modeling'
include { SAVEPARAMS } from '../modules/logging'

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
    pfile_filtered = MAFFILTER(pfile_renamed, params.train)

    //get list of gene names in the input expression table
    //get list of gene names in the QTL table
    //channel with modelable genes based on intersection
    def expr_header = file(params.expr).readLines()[0].tokenize('\t')
    def ensg_cols = expr_header[1..-1]
    expr_ensg = Channel.fromList(ensg_cols)

    qtlindex = GETQTLPICKLE(file(params.qtltable))
    genes = GENELISTPICKLE(qtlindex)
    gene_ch = genes.splitText()
        .map { it.trim() }
        .join(expr_ensg)

    expr_table = PREPROCESSEXPR(pfile_filtered, params.expr, genes)[0]

    // based on modelable gene list create channel of tuples
    //[tissue, ensg, chrom, start, stop]
    params.tis  = file(params.expr).baseName
    gene_info_ch = Channel.fromPath(params.geneInfo)
    
    GETSTARTSTOP(
        gene_info_ch,
        mode_params.upstream,
        mode_params.downstream
    )
    
    gene_info = GETSTARTSTOP.out
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            tuple(row.ENSG, row.CHROM, row.start, row.end)
        }
        .join(gene_ch)
        .map { ensg, chrom, start, stop ->
            [ params.tis, ensg, chrom, start, stop ]
        }

    switch(mode_params.ldmode){
        case "ldnone":
            variant_res = GETVARS(gene_info, pfile_filtered, params.qtltable, qtlindex, false, "", "")
            break
        default:
            variant_res = GETVARS(gene_info, pfile_filtered, params.qtltable, qtlindex, true,  ld_params.ldWindow, ld_params.ldR)
            break
    }

    variants = variant_res
        .map { pgen, psam, pvar ->
            def filename = pgen.baseName  
            def match = filename =~ /^(.+)_(ENSG\d+)$/
            def tis = match[0][1] 
            def ensg = match[0][2]  
            [tis, ensg, pgen, psam, pvar]
        }
            
    models = FITMODEL(variants, expr_table, qtlindex, mode_params.model, mode_params.threshold, params.train)
    SAVEPARAMS(mode_params, ld_params)

}