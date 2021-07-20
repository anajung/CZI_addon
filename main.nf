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
    container 'neherlab/nextclade:latest'
    cpus 1
    memory '1 GB'
    publishDir params.outdir

    input:
    path combined_fa

    output:
    path 'nextclade_lineage.tsv'

    shell:
    //docker run --rm -u 1000 --volume="$PWD:/seq" neherlab/nextclade nextclade --input-fasta '/seq/combined.fa' --output-tsv '/seq/nextclade_lineage.tsv'
    '''
    nextclade --input-fasta !{combined_fa} --output-tsv nextclade_lineage.tsv
    '''
}

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
    //must fix this big thing...

    shell:
    template 'join_lineage.py'

}

process filter_fa {
    container 'quay.io/biocontainers/biopython:1.78'
    cpus 1
    memory '1 GB'
    publishDir params.outdir

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

workflow {
    //vcfdata=channel.fromPath( params.vcf ).map(vcf -> [vcf, vcf.simpleName])
    //vcfConvert(vcfdata)
    //snpEff(vcfConvert.out.vcf_converted)
    //snpSift(snpEff.out.vcf_annotated)

    combinedfadata=channel.fromPath( params.combinedfa )
    //pangolin(combinedfadata)
    //nextClade(combinedfadata)
    //pangolin_lineage=channel.fromPath('/Users/anajung/Documents/HandleyLab/ana_Pangolin_lineage.csv')
    //nextclade_lineage=channel.fromPath('/Users/anajung/Documents/HandleyLab/ana_nextclade_lineage.tsv')
    //python=channel.fromPath('/Users/anajung/Documents/HandleyLab/scripts/CZI_addon/join_lineage.py')
    //joinLineage(pangolin_lineage, nextclade_lineage)
    filter_fa(combinedfadata)
    augur(filter_fa.out)

}