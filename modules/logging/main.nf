import groovy.json.JsonOutput

process SAVEPARAMS {
    publishDir params.outdir, mode: 'copy'
    cache false
    
    input:
    val mode_params
    val ld_params
    
    output:
    path "parameters.json"
    
    script:
    def exclude_keys = ['modes', 'ldstrict', 'ldmed', 'ldlax','ldnone'] // Add keys you want to exclude
    def filtered_params = params.findAll { key, value -> !(key in exclude_keys) }
    def output_data = filtered_params + mode_params + ld_params
    def params_json = JsonOutput.prettyPrint(JsonOutput.toJson(output_data))
    """
    echo '${params_json}' > parameters.json
    """
}