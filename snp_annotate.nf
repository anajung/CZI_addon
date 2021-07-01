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

// process other {
//     shell:
        // if(field_type == 'CZ')
        //     //filed_list = 'POS REF ALT DP4 "ANN[0].EFFECT" "ANN[0].GENE" "ANN[0].HGVS_P"'
        //     '''
        //     java -jar /usr/local/share/snpsift-4.3.1t-3/SnpSift.jar extractFields !{ann_vcf} POS REF ALT DP4 "ANN[0].EFFECT" "ANN[0].GENE" "ANN[0].HGVS_P" > !{vcfID}.snpsift.txt
        //     '''
//     else if(field_type == 'DeepVariant')
//         //filed_list = 'POS REF ALT "ANN[0].EFFECT" "ANN[0].GENE" "GEN[0].VAF" "ANN[0].HGVS_P"'
//         '''
//         java -jar /usr/local/share/snpsift-4.3.1t-3/SnpSift.jar extractFields !{ann_vcf} POS REF ALT "ANN[0].EFFECT" "ANN[0].GENE" "GEN[0].VAF" "ANN[0].HGVS_P" > !{vcfID}.snpsift.txt
//         '''
//     else if(field_type == 'LoFreq')
//         //filed_list = 'POS REF ALT AF "ANN[0].EFFECT" "ANN[0].GENE" "ANN[0].HGVS_P"'
//         '''
//         java -jar /usr/local/share/snpsift-4.3.1t-3/SnpSift.jar extractFields !{ann_vcf} POS REF ALT AF "ANN[0].EFFECT" "ANN[0].GENE" "ANN[0].HGVS_P" > !{vcfID}.snpsift.txt
//         '''
//     //else
//     //    field_list = params.fields

// }

workflow {
    vcfdata=channel.fromPath( params.input ).map(vcf -> [vcf, vcf.simpleName])
    vcfConvert(vcfdata)
    snpEff(vcfConvert.out.vcf_converted)
    snpSift(snpEff.out.vcf_annotated)
}