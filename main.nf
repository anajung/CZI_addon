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

process vcfConvert {
    container 'anajung/pandas'
    cpus 1
    memory '1 GB'
    publishDir params.outdir+'/snpEff_results', mode: 'copy'

    input:
    tuple path(vcf), val(vcfID)
    //vcfID will be retained from input to output 

    output:
    tuple path('*.converted.vcf'), val(vcfID), emit: vcf_converted

    shell:
    '''
    gunzip -f !{vcf}
    sed 's/MN908947.3/NC_045512.2/g' !{vcfID}.vcf > !{vcfID}.converted.vcf
    '''
}

process snpEff {
    container 'quay.io/biocontainers/snpeff:5.0--hdfd78af_1'
    cpus 1
    memory '1 GB'
    publishDir params.outdir+'/snpEff_results', mode: 'copy'

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
    publishDir params.outdir, mode: 'copy'

    input:
    tuple path(ann_vcf), val(vcfID)

    output:
    path '*'

    shell:
    '''
    java -jar /usr/local/share/snpsift-4.3.1t-3/SnpSift.jar extractFields !{ann_vcf} POS REF ALT DP4 "ANN[0].EFFECT" "ANN[0].GENE" "ANN[0].HGVS_P" > !{vcfID}.snpsift.txt
    '''   
    
}

process add_fasta {
    container 'anajung/pandas'
    cpus 1
    publishDir params.outdir, mode: 'copy'

    input:
    path combined_fa
    path references_fa

    output:
    path 'full_combined.fa'

    shell:
    '''
    cat !{combined_fa} !{references_fa} > full_combined.fa
    '''

}
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
    container 'nextstrain/nextclade'
    cpus 4
    memory '6 GB'
    publishDir params.outdir, mode: 'copy'

    input:
    path combined_fa
    path nextclade_root_seq
    path nextclade_tree_json
    path nextclade_qc
    path nextclade_genemap

    output:
    path '*nextclade_lineage.tsv'

    shell:
    '''
    nextclade --input-fasta !{combined_fa} --input-root-seq NC_045512.2.fasta --input-tree tree.json --input-qc-config qc.json --input-gene-map genemap.gff --output-tsv nextclade_lineage.tsv
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
    vcfdata=channel.fromPath( params.vcf ).map(vcf -> [vcf, vcf.simpleName])
    vcfConvert(vcfdata)
    snpEff(vcfConvert.out.vcf_converted)
    snpSift(snpEff.out.vcf_annotated)

    combinedfadata=channel.fromPath( params.combinedfa ).collect()
    tree_references_fasta = channel.fromPath( params.references_fasta ).collect()
    add_fasta(combinedfadata, tree_references_fasta)
    pangolin(add_fasta.out)
    //nextclade_data=channel.fromPath('/Users/anajung/Documents/HandleyLab/scripts/CZI_addon/nextclade_data/*').collect()
    nextclade_root_seq = channel.fromPath( params.root_seq ).collect()
    nextclade_tree_json = channel.fromPath( params.tree_json ).collect()
    nextclade_qc = channel.fromPath( params.qc_json ).collect()
    nextclade_genemap = channel.fromPath( params.genemap ).collect()
    nextClade(add_fasta.out, nextclade_root_seq, nextclade_tree_json, nextclade_qc, nextclade_genemap)
    //nextcladedata=channel.fromPath( params.nextclade )
    joinLineage(pangolin.out, nextClade.out)
    filter_fa(add_fasta.out)
    augur(filter_fa.out)
}