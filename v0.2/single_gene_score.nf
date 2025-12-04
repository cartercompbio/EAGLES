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


params.predixcan = [
    upstream: 1000000,
    downstream: 1000000,
    threshold: 1,
    model: "elasticnet",
    ldmode: "ldnone"
]

params.predixcanLDStrict = [
    upstream: 1000000,
    downstream: 1000000,
    threshold: 1,
    model: "elasticnet",
    ldmode: "ldstrict"
]

params.randomforest = [
    upstream: 1000000,
    downstream: 1000000,
    threshold: 1,
    model: "rf",
    ldmode: "ldlax"
]

params.flipAllele = [
    upstream: 1000000,
    downstream: 1000000,
    threshold: 1,
    model: "flipallele",
    ldmode: "ldstrict"
]

params.magma = [
    upstream: 5000,
    downstream: 1500,
    threshold: 0.999,
    model: "pcr",
    ldmode: "ldnone"
]

params.magmaLDStrict = [
    upstream: 5000,
    downstream: 1500,
    threshold: 0.999,
    model: "pcr",
    ldmode: "ldstrict"
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
    switch(params.mode){
        case "predixcan":
            mode_params=params.predixcan
            break
        case "magma":
            mode_params=params.magma
            break
        case "randomforest":
            mode_params=params.randomforest
            break
        case "flipAllele":
            mode_params = params.flipAllele
            break
        case "predixcanLDStrict":
            mode_params = params.predixcanLDStrict
            break
        case "magmaLDStrict":
            mode_params = params.magmaLDStrict
            break
        default:
            def available_modes = ['predixcan', 'magma', 'randomforest', 'flipAllele', 'predixcanLDStrict', 'magmaLDStrict']
            error "Unknown mode: '${params.mode}'. Available modes: ${available_modes.join(', ')}"
    }
    
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
        .take(100) // Limit to 100 genes for now
    
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
        .collate(10)
        .take(5)        
        
    switch(mode_params.ldmode){
        case "ldnone":
            var_results = GETVARS_BATCH(filtered_tis_gene_ch, pfile_renamed)
            break
        default:
            var_results = GETLDFILTERVARS_BATCH(filtered_tis_gene_ch, pfile_renamed, ld_params.ldWindow, ld_params.ldR)  
            break
    }
    variants = var_results
        .map { pgens, pvars, psams ->
            def pairs = []
            pgens.each { pgen ->
                // Extract basename (e.g., "Brain_ENSG001" from "Brain_ENSG001.pgen")
                def basename = pgen.baseName

                // Extract tis and ensg from basename
                // Assuming format: tissue_ENSG... where ENSG marks the start of gene ID
                def ensgMatch = (basename =~ /^(.+?)_(ENSG\d+.*)$/)
                if (ensgMatch) {
                    def tis = ensgMatch[0][1]
                    def ensg = ensgMatch[0][2]

                    // Find matching pvar and psam files using basename
                    def pvar = pvars.find { it.baseName == basename }
                    def psam = psams.find { it.baseName == basename }

                    if (pvar && psam) {
                        pairs << tuple(tis, ensg, pgen, psam, pvar)
                    }
                }
            }
            return pairs
        }
        .flatMap() //[tis, ensg, pgen, psam, pvar]
        
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
    python ${projectDir}/bin/feature.py --pgen ${pgen} --psam ${psam} --pvar ${pvar} --method ${mode} --thres ${thres} --output "${tis}_${ensg}"
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

process GETVARFEATURES_BATCH {
    cpus 4
    memory 32.GB
    publishDir params.outdir + '/features'
    errorStrategy 'ignore'
    
    input:
    val(batch)  // List of tuples: [(tis, ensg, chrom, start, stop), ...]
    tuple path(pgen), path(pvar), path(psam)
    val(mode)
    val(thres)
    
    output:
    tuple path("*_feats.tsv"), path("*_loadings.tsv"), optional: true
    
    script:
    """
    #!/bin/bash
    
    # Process each item in the batch
    cat << 'EOF' > batch_items.txt
${batch.collect { tis, ensg, chrom, start, stop -> "${tis}\t${ensg}\t${chrom}\t${start}\t${stop}" }.join('\n')}
EOF

    while IFS=\$'\\t' read -r tis ensg chrom start stop; do
        echo "Processing \${tis} - \${ensg}"
        
        python ${projectDir}/bin/qtl_filter.py \\
            --qtl-folder ${params.gtexQTLfolder} \\
            --index-folder ${params.gtexQTLindexFolder} \\
            --gene "\${ensg}" \\
            --tis "\${tis}" \\
            --output "temp_\${tis}_\${ensg}.txt"
        
        if [ -s "temp_\${tis}_\${ensg}.txt" ]; then
            plink2 \\
                --pfile ${pgen.baseName} \\
                --chr "\${chrom}" \\
                --from-bp "\${start}" \\
                --to-bp "\${stop}" \\
                --extract "temp_\${tis}_\${ensg}.txt" \\
                --force-intersect \\
                --make-pgen \\
                --out "\${tis}_\${ensg}"
                
            if [ -s "\${tis}_\${ensg}.pvar" ]; then
                python ${projectDir}/bin/feature.py \\
                    --pgen "\${tis}_\${ensg}.pgen" \\
                    --psam "\${tis}_\${ensg}.psam" \\
                    --pvar "\${tis}_\${ensg}.pvar" \\
                    --method ${mode} \\
                    --thres ${thres} \\
                    --output "\${tis}_\${ensg}"
            else
                echo "No variants for \${tis}_\${ensg}"
            fi
        else
            echo "temp_\${tis}_\${ensg}.txt does not exist or is empty"
        fi
        
        # Clean up intermediate files
        rm -f temp_\${tis}_\${ensg}.txt
        rm -f \${tis}_\${ensg}.pgen \${tis}_\${ensg}.pvar \${tis}_\${ensg}.psam
        
    done < batch_items.txt
    """
}

process GETVARS_BATCH {
    cpus 2
    memory 32.GB
    publishDir params.outdir + '/vars'
    //errorStrategy 'ignore'
    
    input:
    val(batch)  // batch is a list of tuples
    tuple path(pgen), path(pvar), path(psam)
    
    output:
    tuple path("*.pgen"), path("*.pvar"), path("*.psam"), optional: true
    
    script:
    """
    # Process each tuple in the batch
    for item in ${batch.collect { "'${it[0]}|${it[1]}|${it[2]}|${it[3]}|${it[4]}'" }.join(' ')}; do
        IFS='|' read -r tis ensg chrom start stop <<< "\$item"
        
        python ${projectDir}/bin/qtl_filter.py \
            --qtl-folder ${params.gtexQTLfolder} \
            --index-folder ${params.gtexQTLindexFolder} \
            --gene "\$ensg" \
            --tis "\$tis" \
            --output "temp_\${tis}_\${ensg}.txt"
        
        if [ -s "temp_\${tis}_\${ensg}.txt" ]; then
            plink2 \
                --pfile ${pgen.baseName} \
                --chr "\$chrom" \
                --from-bp "\$start" \
                --to-bp "\$stop" \
                --extract "temp_\${tis}_\${ensg}.txt" \
                --force-intersect \
                --make-pgen \
                --out "\${tis}_\${ensg}"
        else
            echo "temp_\${tis}_\${ensg}.txt does not exist or is empty"
        fi
    done
    """
}

//TODO: errors can still arise because some of the gtex variants are missing from the tcga pfiles being used as ld references
//      maybe updating to the tcga wgs would help?

process GETLDFILTERVARS_BATCH {
    cpus 2
    memory 32.GB
    publishDir params.outdir + '/vars'
    //errorStrategy 'ignore'
    
    input:
    val(batch)  // batch is a list of tuples
    tuple path(pgen), path(pvar), path(psam)
    val(window)
    val(r2)
    
    output:
    tuple path("*.pgen"), path("*.pvar"), path("*.psam"), optional: true
    
    script:
    """
    # Process each tuple in the batch
    for item in ${batch.collect { "'${it[0]}|${it[1]}|${it[2]}|${it[3]}|${it[4]}'" }.join(' ')}; do
        IFS='|' read -r tis ensg chrom start stop <<< "\$item"
        
        plink2 \
            --pfile ${pgen.baseName} \
            --chr "\$chrom" \
            --from-bp "\$start" \
            --to-bp "\$stop" \
            --make-pgen \
            --out "\${tis}_\${ensg}_temp"
        
        python ${projectDir}/bin/qtl_filter.py \
            --qtl-folder ${params.gtexQTLfolder} \
            --index-folder ${params.gtexQTLindexFolder} \
            --gene "\$ensg" \
            --tis "\$tis" \
            --pvar "\${tis}_\${ensg}_temp.pvar" \
            --output "temp_\${tis}_\${ensg}.txt"            
        
            
        if [ -s "temp_\${tis}_\${ensg}.txt" ]; then
            num_lines=\$(wc -l < "temp_\${tis}_\${ensg}.txt" )
            
            # exactly 1 eqtl
            if [ "\$num_lines" -eq 1 ]; then
                plink2 \
                    --pfile ${pgen.baseName} \
                    --chr "\$chrom" \
                    --from-bp "\$start" \
                    --to-bp "\$stop" \
                    --extract "temp_\${tis}_\${ensg}.txt" \
                    --force-intersect \
                    --make-pgen \
                    --out "\${tis}_\${ensg}"  
                    
            # more than 1 eqtl, run ld prune first
            else
                plink2 \
                    --indep-pairwise ${window} ${r2} \
                    --pfile ${params.europfile}  \
                    --extract "temp_\${tis}_\${ensg}.txt" \
                    --out "\${tis}_\${ensg}"
                    
                plink2 \
                    --pfile ${pgen.baseName} \
                    --chr "\$chrom" \
                    --from-bp "\$start" \
                    --to-bp "\$stop" \
                    --extract "\${tis}_\${ensg}.prune.in" \
                    --force-intersect \
                    --make-pgen \
                    --out "\${tis}_\${ensg}"                

            fi            
            
        else
             echo "no eqtls found for \${tis} and \${ensg}"
        fi
        
        rm -f *_temp.*
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
    tuple val(tis), val(ensg), path(features), path(expression), path(covariates), path(qtl)
    val(model_type)
    val(thres)

    output:
    tuple val(tis), val(ensg), path("*.pkl"), path("*_scores.tsv")
    //tuple val(tis), val(ensg), path("*.pkl"), optional: true  // in order to still save model

    script:
    """
    python ${projectDir}/bin/fit_and_score.py \
        --features ${features} \
        --expression ${expression} \
        --covariates ${covariates} \
        --model ${model_type} \
        --gene ${ensg} \
        --qtl ${qtl} \
        --samples ${params.train} \
        --thres ${thres} \
        --out-prefix "${tis}_${ensg}"
    """
}
