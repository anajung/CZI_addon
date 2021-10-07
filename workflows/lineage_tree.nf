#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { pangolin; nextClade; joinLineage } from '../modules/lineage'
include { filter_fa; augur} from '../modules/tree'

workflow lineage {
    take:
        assembly
    main:
        pangolin(assembly)
        nextClade(assembly)
        joinLineage(pangolin.out, nextClade.out)
}

workflow tree {
    take:
        assembly
    main:
        filter_fa(combinedfadata)
        augur(filter_fa.out)
}

workflow lineage_tree{
    assembly = channel.fromPath( params.combinedfa ).collect()
    lineage(assembly)
    tree(assembly)
}