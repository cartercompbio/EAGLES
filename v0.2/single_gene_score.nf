// GENERAL PARAMETERS
//params.geneInfo = "/cellar/shared/carterlab/projects/eagle/grch38_gene_tss.tsv"
params.geneInfo = "/cellar/users/domeyer/EAGLE/test_expr/clean_gene_info.tsv"
//params.europfile = "/carter/controlled/dbGaP/phs000178_TCGA/TOPMED_TCGA/plink2_eur_only/tcga.common.european.noimmunecancers"
params.europfile = "/cellar/users/domeyer/EAGLE/test_expr/ld_reference/GTEx.qc_passed.EUR"
params.gtexQTLfolder = "/cellar/users/domeyer/data/gtex/cis_eqtls/GTEx_EUR_slope_tables_by_ENSG_rsid_snp/by_tissue_type"
params.gtexQTLindexFolder = "/cellar/users/domeyer/data/gtex/cis_eqtls/GTEx_EUR_slope_tables_by_ENSG_rsid_snp/by_tissue_type_index"
//params.heritability = "/cellar/users/domeyer/EAGLE/test_expr/tissue_gene_heritability_0_01.tsv"
params.heritability = "/cellar/users/domeyer/EAGLE/test_expr/tissue_gene_heritability_no_predixcan_missing_snps.tsv"

// TCGA PARAMETERS
//params.pfile = "/carter/controlled/dbGaP/phs000178_TCGA/TOPMED_TCGA/plink2/tcga.common"
//params.expressionfolder = "/cellar/users/domeyer/data/tcga/expr_cn_by_ensg/expression"
//params.covariates  = "/cellar/users/nopopko/projects/eagles/covariates_test.csv"

// GTEX PARAMETERS
params.pfile = "/cellar/users/nopopko/projects/eagles/GTEx_plinkqc/qc_output/GTEx.qc_passed"
params.expressionfolder = "/cellar/users/domeyer/data/gtex/expression/by_tissue"
params.covariates  = "/cellar/shared/carterlab/projects/eagle/v0.2/gtex_covar/gtex.eigenvec"
params.train = "/cellar/users/domeyer/EAGLE/test_expr/eur_train_ids.txt"

params.mode = "predixcan" //default MODE

// keep track of different runs and avoid overwriting
def timestamp = new Date().format('MMM-dd-yyyy-HH.mm')
params.outdir = "/cellar/shared/carterlab/projects/eagle/v0.2/test_out/${params.mode}_${timestamp}"


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
    
    randomforest: [
        upstream: 1000000,
        downstream: 1000000,
        threshold: 1,
        model: "rf",
        ldmode: "ldlax"
    ],
    
    flipAllele: [
        upstream: 1000000,
        downstream: 1000000,
        threshold: 1,
        model: "flipallele",
        ldmode: "ldstrict"
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
    
    magma: [
        upstream: 5000,
        downstream: 1500,
        threshold: 0.999,
        model: "pcr",
        ldmode: "ldnone"
    ],
    
    magmaLDStrict: [
        upstream: 5000,
        downstream: 1500,
        threshold: 0.999,
        model: "pcr",
        ldmode: "ldstrict"
    ]
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
    path(gene_info)
    val upstream
    val downstream
    
    output:
    path("*.tsv")
    
    
    script:
    """
    python ${projectDir}/bin/model_window.py \
        --path ${gene_info} \
        --upstream ${upstream} \
        --downstream ${downstream} \
        --output "gene_info"
    """
}

process GETVARS_BATCH {
    cpus 1
    memory { 8.GB * task.attempt }
    maxRetries 4
    maxForks 50
    publishDir params.outdir + '/vars'
    //errorStrategy 'ignore'
    
    input:
    tuple val(tis_list), val(ensg), val(chrom), val(start), val(stop)
    tuple path(pgen), path(pvar), path(psam)
    val(prune)
    val(window)
    val(r2)
    
    output:
    tuple path("*_${ensg}.pgen"), path("*_${ensg}.psam"), path("*_${ensg}.pvar"), optional: true

    
    script:
    """
    # Process each tissue for this gene/region combination
    #for tis in ${tis_list.collect { "'${it}'" }.join(' ')}; do
    for tis in ${tis_list.join(' ')}; do
        
        plink2 \
            --pfile ${pgen.baseName} \
            --chr "${chrom}" \
            --from-bp "${start}" \
            --to-bp "${stop}" \
            --make-pgen \
            --out "\${tis}_${ensg}_temp"
        
        python ${projectDir}/bin/qtl_filter.py \
            --qtl-folder ${params.gtexQTLfolder} \
            --index-folder ${params.gtexQTLindexFolder} \
            --gene "${ensg}" \
            --tis "\${tis}" \
            --pvar "\${tis}_${ensg}_temp.pvar" \
            --output "temp_\${tis}_${ensg}.txt"            
        
        # any eqtl found
        if [ -s "temp_\${tis}_${ensg}.txt" ]; then
            num_lines=\$(grep -c '' "temp_\${tis}_${ensg}.txt")
            
            # exactly 1 eqtl
            if [ "\$num_lines" -eq 1 ]; then
                plink2 \
                    --pfile ${pgen.baseName} \
                    --chr "${chrom}" \
                    --from-bp "${start}" \
                    --to-bp "${stop}" \
                    --extract "temp_\${tis}_${ensg}.txt" \
                    --force-intersect \
                    --make-pgen \
                    --out "\${tis}_${ensg}"  
                    
            # more than 1 eqtl and ld prune enabled
            elif [ "${prune}" = "true" ]; then
                plink2 \
                    --indep-pairwise ${window} ${r2} \
                    --pfile ${params.europfile}  \
                    --extract "temp_\${tis}_${ensg}.txt" \
                    --out "\${tis}_${ensg}"
                    
                plink2 \
                    --pfile ${pgen.baseName} \
                    --chr "${chrom}" \
                    --from-bp "${start}" \
                    --to-bp "${stop}" \
                    --extract "\${tis}_${ensg}.prune.in" \
                    --force-intersect \
                    --make-pgen \
                    --out "\${tis}_${ensg}"
                    
            # more than 1 eqtl, no ld pruning
            else
                plink2 \
                    --pfile ${pgen.baseName} \
                    --chr "${chrom}" \
                    --from-bp "${start}" \
                    --to-bp "${stop}" \
                    --extract "temp_\${tis}_${ensg}.txt" \
                    --force-intersect \
                    --make-pgen \
                    --out "\${tis}_${ensg}"
            fi        
        #no eqtls found    
        else
             echo "no eqtls found for \${tis} and ${ensg}"
        fi
        
        rm -f *temp*
    done
    
    
    """
}

process FITMODEL{
    cpus 1
    memory 32.GB
    publishDir params.outdir + '/models'
    //errorStrategy 'ignore'
    
    input:
    tuple val(tis), val(ensg), path(pgen), path(psam), path(pvar), path(expression), path(covariates), path(qtl)
    val(model_type)
    val(thres)
    path(train)

    output:
    tuple val(tis), val(ensg), path("*.pkl")
    
    script:
    """
    python ${projectDir}/bin/fit_model.py \\
        --pgen ${pgen} \\
        --psam ${psam} \\
        --pvar ${pvar} \\
        --expression ${expression} \\
        --covariates ${covariates} \\
        --model ${model_type} \\
        --gene ${ensg} \\
        --qtl ${qtl} \\
        --samples ${train} \\
        --thres ${thres} \\
        --output "${tis}_${ensg}.pkl"
    """ 
}

process MODELSCORE {
    cpus 1
    memory 32.GB
    publishDir params.outdir + '/scores'
    
    input:
    tuple val(tis), val(ensg), path(pgen), path(psam), path(pvar), path(model), path(covariates)

    output:
    tuple val(ensg), path("*_scores.tsv")
    
    script:
    """
    python ${projectDir}/bin/model_score.py \\
        --pgen ${pgen} \\
        --psam ${psam} \\
        --pvar ${pvar} \\
        --covariates ${covariates} \\
        --model ${model} \\
        --output "${tis}_${ensg}_scores.tsv"
    """
}

process FIT_AND_SCORE {
    cpus 2
    memory 32.GB
    publishDir params.outdir + '/scores', pattern: '*_scores.tsv'
    publishDir params.outdir + '/models', pattern: '*.pkl'

    input:
    tuple val(tis), val(ensg), path(pgen), path(psam), path(pvar), path(expression), path(covariates), path(qtl)
    val(model_type)
    val(thres)

    output:
    tuple val(tis), val(ensg), path("*.pkl"), path("*_scores.tsv")
    //tuple val(tis), val(ensg), path("*.pkl"), optional: true  // in order to still save model

    script:
    """
    python ${projectDir}/bin/fit_and_score.py \
        --pgen ${pgen} \
        --psam ${psam} \
        --pvar ${pvar} \
        --expression ${expression} \
        --covariates ${covariates} \
        --model ${model_type} \
        --gene ${ensg} \
        --qtl ${qtl} \
        --samples ${params.train} \
        --thres ${thres} \
        --model_output "${tis}_${ensg}.pkl" \
        --score_output "${tis}_${ensg}"
    """
}
