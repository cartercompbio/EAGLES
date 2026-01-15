process remove_chr_prefix {
    publishDir params.qc_dir, mode: 'symlink', overwrite: true

    input:
    tuple val(prefix), path(pgen_pvar_psam)

    output:
    tuple val("${prefix}.noch"), path(outs)

    script:
    outs = [
        "${prefix}.noch.pgen",
        "${prefix}.noch.pvar",
        "${prefix}.noch.psam",
    ]
    """
    cp ${prefix}.pgen ${prefix}.noch.pgen
    cp ${prefix}.psam ${prefix}.noch.psam
    awk 'BEGIN{OFS="\\t"} /^#/{print;next} {sub(/^chr/, "", \$1); print}' ${prefix}.pvar > ${prefix}.noch.pvar
    """
}


process recode_study_var_ids {
    publishDir params.qc_dir, mode: 'symlink', overwrite: true

    input:
    tuple val(prefix), path(pgen_pvar_psam)

    output:
    tuple val("${prefix}.recode"), path(outs)

    script:
    outs = [
        "${prefix}.recode.pgen",
        "${prefix}.recode.pvar",
        "${prefix}.recode.psam",
    ]
    """
    plink2 --pfile ${prefix} \\
        --set-all-var-ids @:#:\\\$r:\\\$a \\
        --new-id-max-allele-len 30 missing \\
        --make-pgen \\
        --out ${prefix}.recode
    """
}

process plink_stats {
    publishDir params.qc_dir, mode: 'symlink', overwrite: true

    input:
    tuple val(prefix), path(pgen_pvar_psam)

    output:
    path "${prefix}.smiss"
    path "${prefix}.vmiss"
    path "${prefix}.het"
    path "${prefix}.acount"
    path "${prefix}.hardy"

    script:
    """
    plink2 --pfile $prefix \\
    --missing \\
    --het \\
    --freq counts \\
    --hardy midp log10 \\
    --out $prefix
    """
}

process ac_gt_snps {
    input:
    tuple val(prefix), path(pgen_pvar_psam)

    output:
    tuple val(prefix), path(outs)

    script:
    outs = [
        "${prefix}.no_ac_gt_snps.pgen",
        "${prefix}.no_ac_gt_snps.pvar",
        "${prefix}.no_ac_gt_snps.psam",
    ]
    """
    awk 'BEGIN {OFS="\t"} !/^#/ && (\$4\$5 == "GC" || \$4\$5 == "CG" || \$4\$5 == "AT" || \$4\$5 == "TA") {print \$3}' \\
    ${prefix}.pvar > ${prefix}.ac_gt_snps

    plink2 --pfile $prefix \\
    --exclude ${prefix}.ac_gt_snps \\
    --make-pgen \\
    --out ${prefix}.no_ac_gt_snps
    """
}

process prune {
    publishDir params.qc_dir, mode: 'symlink', overwrite: true, pattern: '*.{pgen,pvar,psam}'

    input:
    tuple val(prefix), path(pgen_pvar_psam)
    path highld

    output:
    tuple val(prefix), path(outs)
    path "${filt_prefix}.prune.in"

    script:
    filt_prefix = "${prefix}.no_ac_gt_snps"
    outs = [
        "${prefix}.pruned.pgen",
        "${prefix}.pruned.pvar",
        "${prefix}.pruned.psam",
    ]
    """
    plink2 --pfile $filt_prefix \\
    --exclude range $highld \\
    --set-all-var-ids @:#:\\\$r:\\\$a \\
    --rm-dup exclude-all \\
    --new-id-max-allele-len 1000 missing \\
    --indep-pairwise 50 5 0.7 \\
    --out $filt_prefix

    plink2 --pfile $filt_prefix \\
    --extract ${filt_prefix}.prune.in \\
    --make-pgen \\
    --out ${prefix}.pruned
    """
}

process identity_by_descent {
    publishDir params.qc_dir, mode: 'symlink', overwrite: true

    input:
    tuple val(prefix), path(pgen_pvar_psam)

    output:
    path "${prefix}.pruned.king.bin"

    script:
    """
    plink2 --pfile ${prefix}.pruned \\
    --make-king triangle bin4 \\
    --out ${prefix}.pruned
    """
}

process filter_reference {
    input:
    tuple val(ref_prefix), path(ref_pgen_pvar_psam)
    path prune_in

    output:
    tuple val(ref_prefix), path(outs)

    script:
    outs = [
        "${ref_prefix}.pruned.pgen",
        "${ref_prefix}.pruned.pvar",
        "${ref_prefix}.pruned.psam",
    ]
    """
    plink2 --pfile $ref_prefix \\
    --extract $prune_in \\
    --make-pgen \\
    --out ${ref_prefix}.pruned
    """
}

process merge_study_ref {
    publishDir params.qc_dir, mode: 'symlink', overwrite: true
    memory '128GB'

    input:
    tuple val(prefix), path(pgen_pvar_psam)
    tuple val(ref_prefix), path(ref_pgen_pvar_psam)

    output:
    tuple val(merged_prefix), path(outs)

    script:
    merged_prefix = "${prefix}.merged.${ref_prefix}"
    outs = [
        "${merged_prefix}.bed",
        "${merged_prefix}.bim",
        "${merged_prefix}.fam",
    ]
    """
    plink2 --pfile ${prefix}.pruned \\
    --extract-col-cond ${ref_prefix}.pruned.pvar 2 3 '#' \\
    --make-bed \\
    --memory ${task.memory.toMega()} \\
    --out ${prefix}.pruned

    plink2 --pfile ${ref_prefix}.pruned \\
    --extract-col-cond ${prefix}.pruned.pvar 2 3 '#' \\
    --make-bed \\
    --memory ${task.memory.toMega()} \\
    --out ${ref_prefix}.pruned

    plink --bfile ${prefix}.pruned \\
    --bmerge ${ref_prefix}.pruned \\
    --make-bed \\
    --memory ${task.memory.toMega()} \\
    --out $merged_prefix
    """
}

process pca {
    publishDir params.qc_dir, mode: 'symlink', overwrite: true
    cpus 8
    memory '64GB'

    input:
    tuple val(prefix), path(bed_bim_bam)

    output:
    path "${prefix}.eigenvec"

    script:
    """
    plink2 --bfile $prefix \\
    --geno 0.01 \\
    --maf 0.01 \\
    --pca \\
    --out $prefix
    """
}

process plot_stats {
    publishDir params.qc_dir, mode: 'copy', overwrite: true
    cpus 16
    memory '256GB'

    input:
    path smiss
    path vmiss
    path het
    path acount
    path hardy

    output:
    path "*.png"
    path "${params.cohort_name}.vpass"
    path "${params.cohort_name}.spass"

    script:
    """
    plot_stats.py \\
    $params.cohort_name \\
    $smiss \\
    $vmiss \\
    $het \\
    $acount \\
    $hardy
    """
}

process plot_ibd {
    publishDir params.qc_dir, mode: 'copy', overwrite: true
    cpus 1
    memory '16GB'

    input:
    path king_bin

    output:
    path "*.png"

    script:
    """
    plot_id_by_descent.py $king_bin
    """
}

process plot_pca {
    publishDir params.qc_dir, mode: 'copy', overwrite: true
    cpus 1
    memory '8GB'

    input:
    path eigenvec

    output:
    path "*.png"

    script:
    """
    plot_pca.py $params.cohort_name $eigenvec
    """
}
