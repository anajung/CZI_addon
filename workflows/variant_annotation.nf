#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { vcfConvert; snpEff; snpSift } from '../modules/snp'

workflow variant_annotation {
    take:
        vcfdata
    main:
        vcfConvert(vcfdata)
        snpEff(vcfConvert.out.vcf_converted)
        snpSift(snpEff.out.vcf_annotated)

}

