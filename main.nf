#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process vcfConvert {
    container 'quay.io/biocontainers/snpeff:5.0--hdfd78af_1'
    cpus 1
    memory '1 GB'
    publishDir params.outdir

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
    publishDir params.outdir

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

//field_type = params.fields

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
    path 'pangolin_lineage.csv'

    //pangolin --usher combined.fa --outfile usher_Pangolin_lineage.csv
    shell:
    '''
    pangolin --usher !{combined_fa} --outfile pangolin_lineage.csv
    '''
}

process nextClade {
    container 'neherlab/nextclade'
    cpus 1
    memory '1 GB'
    publishDir params.outdir

    input:
    path combined_fa

    output:
    path 'nextclade_lineage.csv'

    shell:
    '''
    nextclade --input-fasta !{combined_fa} --output-csv nextclade_lineage.csv
    '''
}

process joinLineage {
    container 'amancevice/pandas'
    cpus 1
    memory '1 GB'
    publishDir params.outdir

    input:
    path pangolin_file
    path nextclade_file
    path python

    output:
    path 'joined_lineage.csv'
    //must fix this big thing...
    shell:
    '''
    python join_lineage.py
    '''

}

workflow {
    vcfdata=channel.fromPath( params.vcf ).map(vcf -> [vcf, vcf.simpleName])
    vcfConvert(vcfdata)
    snpEff(vcfConvert.out.vcf_converted)
    snpSift(snpEff.out.vcf_annotated)

    combinedfadata=channel.fromPath( params.combinedfa )
    pangolin(combinedfadata)
    nextClade(combinedfadata)
    //pangolin_file=channel.fromPath('/Users/anajung/Documents/HandleyLab/ana_Pangolin_lineage.csv')
    //nextclade_file=channel.fromPath('/Users/anajung/Documents/HandleyLab/ana_nextclade_lineage.tsv')
    //python=channel.fromPath('/Users/anajung/Documents/HandleyLab/scripts/CZI_addon/join_lineage.py')
    //joinLineage(python, pangolin_file, nextclade_file)

}