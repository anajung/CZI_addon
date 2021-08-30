#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def helpMessage() {
  log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run anajung/CZI_addon --combinedfa 'combined.fa' --vcf 'sample-variants/*.vcf.gz' --outdir './out'

    Mandatory arguments:
    --combinedfa                Path to combined FASTA, must be in quotes
    --outdir                    Specify path to output directory

    Optional arguments:
    --vcf                       Path to .vcf.gz's if variant annotation workflow is wanted
    """.stripIndent()
}

// if (params.help) workflow {
//     helpMessage()
//     exit 0
// }

process pangolin {
    container 'staphb/pangolin'
    cpus 1
    memory '1 GB'
    publishDir params.outdir, mode: 'copy'

    input:
    path combined_fa

    output:
    path '*_lineage.csv'

    shell:
    '''
    pangolin --usher !{combined_fa} --outfile pangolin_lineage.csv
    '''
}

process nextClade {
    //container 'neherlab/nextclade:0.14.4-stretch'
    container 'anajung/nextclade'
    cpus 1
    memory '1 GB'
    publishDir params.outdir, mode: 'copy'

    input:
    path combined_fa

    output:
    path '*nextclade_lineage.tsv'

    shell:
    '''
    nextclade -i !{combined_fa} -t nextclade_lineage.tsv
    '''
}

process joinLineage {
    container 'anajung/pandas'
    cpus 1
    memory '1 GB'
    publishDir params.outdir, mode: 'copy'

    input:
    path pangolin_lineage
    path nextclade_lineage

    output:
    path 'joined_lineage.csv'
    path 'filtered_joined_lineage.csv'

    shell:
    template 'join_lineage.py'

}

process filter_fa {
    container 'quay.io/biocontainers/biopython:1.78'
    cpus 1
    memory '1 GB'

    input:
    path combinedfadata

    output:
    path '*.fa'

    shell:
    template 'filter_fasta.py'
}

process augur {
    container 'anajung/nextstrain'
    cpus 1
    memory '1 GB'
    publishDir params.outdir, mode: 'copy'

    input:
    path filtered_fa

    output:
    path '*.nwk'

    shell:
    '''
    augur align -s !{filtered_fa}
    augur tree -a alignment.fasta -o tree.nwk
    '''
}

workflow {
    combinedfadata=channel.fromPath( params.combinedfa ).collect()
    pangolin(combinedfadata)
    //nextClade(combinedfadata)
    nextcladedata=channel.fromPath( params.nextclade )
    joinLineage(pangolin.out, nextcladedata)
    filter_fa(combinedfadata)
    augur(filter_fa.out)
}
