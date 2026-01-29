params.pfile = "/cellar/users/nopopko/projects/eagles/grievous_dbs/geuvadis_aligned"
params.covariates = "${params.pfile}.psam"

plink2_prefix = file(params.pfile).getBaseName()
params.outdir = file("${params.modelfolder}/${plink2_prefix}")

include { RENAMEVARIANTS; GATHERGENO } from '../modules/preprocessing'
include { MODELSCORE } from '../modules/modeling'

def validateParams() {

    if (!params.modelfolder) {
        error "ERROR: --modelfolder parameter is required"
    }
    if (!params.pfile) {
        error "ERROR: --pfile parameter is required"
    }
    
    // Validate outdir
    //def outdir = file(params.outdir)
    if (params.outdir.exists() && params.outdir.isDirectory()) {
        def scoresDir = file("${params.outdir}/scores")
        if (scoresDir.exists() && scoresDir.isDirectory()) {
            error "ERROR: Output directory '${params.outdir}' already contains a 'scores' folder. Please use a different output directory or remove the existing 'scores' folder."
        }
    }
    
    // Validate modelfolder exists
    def modelfolder = file(params.modelfolder)
    if (!modelfolder.exists() || !modelfolder.isDirectory()) {
        error "ERROR: Model folder '${params.modelfolder}' does not exist or is not a directory"
    }
    
    
    // Validate pfile prefix
    def pgen = file("${params.pfile}.pgen")
    def psam = file("${params.pfile}.psam")
    def pvar = file("${params.pfile}.pvar")
    
    def missing = []
    if (!pgen.exists()) missing.add(".pgen")
    if (!psam.exists()) missing.add(".psam")
    if (!pvar.exists()) missing.add(".pvar")
    
    if (missing) {
        error "ERROR: Missing PLINK2 files for prefix '${params.pfile}': ${missing.join(', ')}"
    }
}

workflow {

    validateParams()
    
    //if no covariate file specified, uses the given psam file
    def cov_file = params.covariates ?: "${params.pfile}.psam"


    pfile_renamed = RENAMEVARIANTS(params.pfile) // pgen, pvar, psam
    
    model_ch = Channel
        .fromPath("${params.modelfolder}/*.pkl")
        .map { model_file ->
                def filename = model_file.getBaseName()
                def parts = filename.split('_')
                def ensg = parts[-1]
                def tis = parts[0..-2].join('_') 
                [tis, ensg, model_file]
        }

    geno_ch = pfile_renamed
        .combine(model_ch)
        .map { pgen, pvar, psam, tis, ensg, model ->
            def cov = file(params.covariates)
            tuple(tis, ensg, pgen, pvar, psam, model, cov)
        }
        | GATHERGENO  

    //TODO: need to check if any of the model snps are found in the pfile without raising plink2 error
    
    score_ch = MODELSCORE(geno_ch)
   
}