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
