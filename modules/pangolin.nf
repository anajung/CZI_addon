#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process pangolin {
    container 'staphb/pangolin'
    cpus 1
    memory '1 GB'
    publishDir params.outdir

    input:
    path combined_fa

    output:
    path '*_lineage.csv'

    shell:
    '''
    pangolin --usher !{combined_fa} --outfile pangolin_lineage.csv
    '''
}

workflow {
    combinedfadata=channel.fromPath( params.combinedfa )
    pangolin(combinedfadata)
}