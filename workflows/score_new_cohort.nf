
include { GATHERSNPNAMES; RENAMEVARIANTS; SNPSLICE } from '../modules/preprocessing'
include { MULTISCORE } from '../modules/modeling'
workflow {
    model_ch = Channel
        .fromPath("${params.modelfolder}/*.pkl")
        .map { model_file ->
                def filename = model_file.getBaseName()
                def parts = filename.split('_')
                def ensg = parts[-1]
                def tis = parts[0..-2].join('_') 
                [tis, ensg, model_file]
        }

    snp_ch = GATHERSNPNAMES(model_ch.map { tis, ensg, model_file -> model_file })
        .collectFile(newLine: true, name: 'all_snps.txt', storeDir: "${params.outdir}/snplist")

    pfile_renamed = RENAMEVARIANTS(params.pfile)

    pfile_sliced = SNPSLICE(pfile_renamed, snp_ch)

    scores = MULTISCORE(pfile_sliced)

}