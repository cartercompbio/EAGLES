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
    errorStrategy 'ignore'
    
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
    #given window specified by (chrom, start-stop) determine how many variants are present
    #generate slice of pfile with these variants
    count=\$(awk -v col1=${chrom} -v lower=${start} -v upper=${stop} 'NF >= 2 && \$1 == col1 && \$2 >= lower && \$2 <= upper {c++} END {print c+0}' "${pvar}")
    if [ \$count -eq 0 ]; then
        touch "${ensg}_temp"
    else 
        plink2 \
            --pfile ${pgen.baseName} \
            --chr "${chrom}" \
            --from-bp "${start}" \
            --to-bp "${stop}" \
            --make-pgen \
            --out "${ensg}_temp"
    fi

    for tis in ${tis_list.join(' ')}; do
        #identify tis-associated eqtls from pfile
        python ${projectDir}/bin/qtl_filter.py \
            --qtl-folder ${params.gtexQTLfolder} \
            --index-folder ${params.gtexQTLindexFolder} \
            --gene "${ensg}" \
            --tis "\${tis}" \
            --pvar "${ensg}_temp" \
            --output "temp_\${tis}_${ensg}.txt"            
        
        # any eqtl found
        if [ -s "temp_\${tis}_${ensg}.txt" ]; then
            num_lines=\$(grep -c '' "temp_\${tis}_${ensg}.txt")
            
            # exactly 1 eqtl
            if [ "\$num_lines" -eq 1 ]; then
                plink2 \
                    --pfile "${ensg}_temp" \
                    --extract "temp_\${tis}_${ensg}.txt" \
                    --make-pgen \
                    --out "\${tis}_${ensg}"  
                    
            # more than 1 eqtl and ld prune enabled
            elif [ "${prune}" = "true" ]; then
                plink2 \
                    --indep-pairwise ${window} ${r2} \
                    --pfile ${params.europfile}  \
                    --extract "temp_\${tis}_${ensg}.txt" \
                    --out "\${tis}_${ensg}"
                    
                if [ -s "\${tis}_${ensg}.prune.in" ]; then
                    plink2 \
                        --pfile "${ensg}_temp" \
                        --extract "\${tis}_${ensg}.prune.in" \
                        --make-pgen \
                        --out "\${tis}_${ensg}"
                else
                    echo "no lo-LD eqtls found for \${tis} and ${ensg}"
                fi
                    
            # more than 1 eqtl, no ld pruning
            else
                plink2 \
                    --pfile "${ensg}_temp" \
                    --extract "temp_\${tis}_${ensg}.txt" \
                    --make-pgen \
                    --out "\${tis}_${ensg}"
            fi        
        #no eqtls found    
        else
             echo "no eqtls found for \${tis} and ${ensg}"
        fi
        
    done
    rm -f *temp*
    
    """
}
