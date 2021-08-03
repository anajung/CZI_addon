#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { variant_annotation } './workflows/variant_annotation'
include { lineage; tree } from './workflows/lineage_tree'

workflow lineage_tree {
    lineage(assembly)
    tree(assembly)
}

vcfdata = Channel.fromPath( params.vcf ).map(vcf -> [vcf, vcf.simpleName])
assembly = Channel.fromPath( params.combinedfa ).collect()

workflow {
    variant_annotation(vcfdata)
    lineage_tree(assembly)
}

// workflow {
//     vcfdata=channel.fromPath( params.vcf ).map(vcf -> [vcf, vcf.simpleName])
//     vcfConvert(vcfdata)
//     snpEff(vcfConvert.out.vcf_converted)
//     snpSift(snpEff.out.vcf_annotated)
//     combinedfadata=channel.fromPath( params.combinedfa )
//     pangolin(combinedfadata)
//     nextClade(combinedfadata)
//     joinLineage(pangolin.out, nextClade.out)
//     filter_fa(combinedfadata)
//     augur(filter_fa.out)

// }