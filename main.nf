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

if (params.help) workflow {
    helpMessage()
    exit 0
}

process vcfConvert {
    container 'anajung/pandas'
    cpus 1
    memory '1 GB'
    publishDir params.outdir+'/snpEff_results'

    input:
    tuple path(vcf), val(vcfID)
    //vcfID will be retained from input to output 

    output:
    tuple path('*.converted.vcf'), val(vcfID), emit: vcf_converted

    shell:
    '''
    gunzip !{vcf}
    sed 's/MN908947.3/NC_045512.2/g' !{vcfID}.vcf > !{vcfID}.converted.vcf
    '''
}

process snpEff {
    container 'quay.io/biocontainers/snpeff:5.0--hdfd78af_1'
    cpus 1
    memory '1 GB'
    publishDir params.outdir+'/snpEff_results'

    input:
    tuple path(converted_vcf), val(vcfID)
    
    output:
    tuple path('*.ann.vcf'), val(vcfID), emit: vcf_annotated
    path '*.txt'
    path '*.html'

    shell:
    '''
    java -Xmx8g -jar /usr/local/share/snpeff-5.0-1/snpEff.jar NC_045512.2 !{vcfID}.converted.vcf > !{vcfID}.ann.vcf
    mv snpEff_genes.txt !{vcfID}_snpEff_genes.txt
    mv snpEff_summary.html !{vcfID}_snpEff_summary.html
    '''
}

process snpSift {
    container 'quay.io/biocontainers/snpsift:4.3.1t--hdfd78af_3'
    cpus 1
    memory '1 GB'
    publishDir params.outdir

    input:
    tuple path(ann_vcf), val(vcfID)

    output:
    path '*'

    shell:
    '''
    java -jar /usr/local/share/snpsift-4.3.1t-3/SnpSift.jar extractFields !{ann_vcf} POS REF ALT DP4 "ANN[0].EFFECT" "ANN[0].GENE" "ANN[0].HGVS_P" > !{vcfID}.snpsift.txt
    '''   
    
}

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

process nextClade {
    container 'neherlab/nextclade:0.14.4-stretch'
    cpus 1
    memory '1 GB'
    publishDir params.outdir

    input:
    path combined_fa

    output:
    path '*nextclade_lineage.tsv'

    shell:
    '''
    nextclade --input-fasta !{combined_fa} --output-tsv nextclade_lineage.tsv
    '''
}

process joinLineage {
    container 'anajung/pandas'
    cpus 1
    memory '1 GB'
    publishDir params.outdir

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
    publishDir params.outdir

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

vcfdata=channel.fromPath( params.vcf ).map(vcf -> [vcf, vcf.simpleName])
combinedfadata=channel.fromPath( params.combinedfa ).collect()


workflow {
    if ( params.vcf )
        vcfConvert(vcfdata)
        snpEff(vcfConvert.out.vcf_converted)
        snpSift(snpEff.out.vcf_annotated)
    pangolin(combinedfadata)
    nextClade(combinedfadata)
    joinLineage(pangolin.out, nextClade.out)
    filter_fa(combinedfadata)
    augur(filter_fa.out)
}
