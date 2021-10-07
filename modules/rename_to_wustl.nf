#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process mgi_to_wustl {
    container 'anajung/pandas'
    cpus 1
    memory '1 GB'
    publishDir params.outdir

    input:
    path samplemap

    output:
    path 'Samplemap.csv'
    path 'rename_guide.txt', emit: rename_guide

    shell:
    template 'samplemap_to_wustl.py'
    
}

process convertName {
    container 'anajung/pandas'
    cpus 1
    memory '1 GB'
    publishDir params.outdir

    input:
    path rename_guide
    path fastq

    output:
    path '*.fastq.gz'

    shell:
    '''
    eval "$(sed 's/^/mv /g' rename_guide.txt )"
    '''
    
}
workflow {
    samplemap = channel.fromPath( params.sample )
    fastq = channel.fromPath( params.fastq )
    mgi_to_wustl(samplemap)
    convertName(mgi_to_wustl.out.rename_guide, fastq)

}