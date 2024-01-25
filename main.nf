#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { VEP_ANNOTATE } from "${baseDir}/modules/VEP/main"
include { VCF2MAF } from "${baseDir}/modules/vcf2maf/main"
include { MAFTOOLS } from "${baseDir}/modules/maftools/main"
include { ANNOVAR_ANNOTATE } from "${baseDir}/modules/annovar/main"
include { ANNOVAR2MAF } from "${baseDir}/nextflow_modules/modules/annovar2maf/main"
include { SEQUENZA_CNAqc } from "${baseDir}/modules/Sequenza_CNAqc/main"
//include { BCFTOOLS_SPLIT_VEP } from "${baseDir}/modules/bcftools/main"
//include { JOIN_TABLES } from "${baseDir}/modules/join_tables/main"
//include { SEQUENZA_CNAqc } from "${baseDir}/modules/Sequenza_CNAqc/main"
//include { SUBCLONAL_DECONVOLUTION } from "${baseDir}/subworkflows/subclonal_deconvolution/main"


workflow {

  input_vcf = Channel.fromPath(params.samples).
      splitCsv(header: true).
      map{row ->
        tuple(row.dataset.toString(), row.patient.toString(), row.sample.toString(), file(row.vcf))}
  

  input_sequenza = Channel.fromPath(params.samples).
    splitCsv(header: true).
    map{row ->
     tuple(row.patient.toString(), row.sample.toString(), row.sex.toString(), file(row.seqz), file(row.vcf))}

vep_out = VEP_ANNOTATE(input_vcf) 
annovar_out = ANNOVAR_ANNOTATE(input_vcf)
vcf2maf_out = VCF2MAF(vep_out)
maf_out = MAFTOOLS(vcf2maf_out.groupTuple(by: 0))
annovar2maf_out = ANNOVAR2MAF(annovar_out)
//BCFTOOLS_SPLIT_VEP(vep_output)
//PLATYPUS_CALL_VARIANTS(input_multisample.groupTuple(by: [0,3,4]))
//JOINT_TABLE()
SEQUENZA_CNAqc(input_sequenza)
//SUBCLONAL_DECONVOLUTION(joint_table)
}
