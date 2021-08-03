#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process joinLineage {
    container 'anabugseq/pandas:latest'
    cpus 1
    memory '1 GB'
    publishDir params.outdir

    input:
    path pangolin_lineage
    path nextclade_lineage

    output:
    path 'joined_lineage.csv'
    path 'filtered_joined_lineage.csv'
    //must fix this big thing...

    shell:
    template 'join_lineage.py'

}

workflow {
    pangolin_lineage = channel.fromPath( params.pangolin )
    nextclade_lineage = channel.fromPath( params.nextclade )
    joinLineage(pangolin_lineage, nextclade_lineage)
}