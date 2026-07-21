process FITMODEL{
    publishDir params.outdir + '/models'
    label 'EAGLES_MODEL'
    
    input:
    tuple val(tis), val(ensg), path(pgen), path(psam), path(pvar), path(expression), path(covariates, stageAs: 'covariates*'), path(qtl), path(qtl_index)
    val(model_type)
    val(thres)
    path(train)

    output:
    tuple val(tis), val(ensg), path("*.pkl"), optional: true
    
    script:
    def cov_arg = covariates instanceof List && covariates.isEmpty() ? "" : "--covariates ${covariates}"
    """
    variant_count=\$(grep -c '^[^#]' ${pvar})
    
    if [ "\$variant_count" -gt 0 ]; then
        python ${projectDir}/bin/fit_model.py \\
            --pgen ${pgen} \\
            --psam ${psam} \\
            --pvar ${pvar} \\
            --expression ${expression} \\
            ${cov_arg} \\
            --model ${model_type} \\
            --gene ${ensg} \\
            --qtl ${qtl} \\
            --qtl-index ${qtl_index} \\
            --samples ${train} \\
            --thres ${thres} \\
            --output "${tis}_${ensg}.pkl"
            
        # output is only optional for non-pcr model types
        if [ "${model_type}" != "pcr" ] && [ ! -f "${tis}_${ensg}.pkl" ]; then
            echo "Error: Expected output file ${tis}_${ensg}.pkl not created"
            exit 1
        fi
        
    else
        echo "Skipping ${tis}_${ensg}: Only \$variant_count variant(s) found in ${pvar}"
    fi
    """ 
}

process MODELSCORE {
    publishDir params.outdir + '/scores'
    label 'EAGLES_MODEL'
    
    input:
    tuple val(tis), val(ensg), path(pgen), path(psam), path(pvar), path(model), path(covariates, stageAs: 'covariates*')

    output:
    tuple val(tis), val(ensg), path("*_scores.tsv")
    
    script:
    def cov_arg = covariates instanceof List && covariates.isEmpty() ? "" : "--covariates ${covariates}"

    """
    python ${projectDir}/bin/model_score.py \\
        --pgen ${pgen} \\
        --psam ${psam} \\
        --pvar ${pvar} \\
        ${cov_arg} \\
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
    path(train)

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
        --samples ${train} \
        --thres ${thres} \
        --model_output "${tis}_${ensg}.pkl" \
        --score_output "${tis}_${ensg}"
    """
}
