// Instructions to run:
// for chr in {1..22}; do
//     nextflow run plinkqc.nf -c configs/chr${chr}.config
// done

params.highld = "../data/high-LD-regions-hg38-GRCh38.bed"
// copied from https://github.com/meyer-lab-cshl/plinkQC/blob/master/inst/extdata/high-LD-regions-hg38-GRCh38.bed

include { remove_chr_prefix; recode_study_var_ids; plink_stats; ac_gt_snps; prune; identity_by_descent; filter_reference; merge_study_ref; pca; plot_stats; plot_ibd; plot_pca } from '../modules/plinkqc'

workflow {
    study = Channel.fromFilePairs("${params.study}.{pgen,pvar,psam}", size: 3)
    ref = Channel.fromFilePairs("${params.ref}.{pgen,pvar,psam}", size: 3)
    
    noch_study = remove_chr_prefix(study)
    
    recoded_study = recode_study_var_ids(noch_study)

    stats = plink_stats(recoded_study)
    plot_stats(stats)
    
    no_ac_gt = ac_gt_snps(recoded_study)
    (pruned, prune_in) = prune(no_ac_gt, params.highld)
    ibd = identity_by_descent(pruned)
    plot_ibd(ibd)

    filter_ref = filter_reference(ref, prune_in)
    merged = merge_study_ref(pruned, filter_ref)
    pca = pca(merged)
    plot_pca(pca)

}
