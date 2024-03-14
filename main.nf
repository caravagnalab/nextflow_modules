#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { VARIANT_ANNOTATION } from "${baseDir}/subworkflows/variant_annotation/main"
include { SUBCLONAL_DECONVOLUTION } from "${baseDir}/subworkflows/subclonal_deconvolution/main"
include { QC } from "${baseDir}/subworkflows/CNAqc/main"
include { PLATYPUS_CALL_VARIANTS } from "${baseDir}/modules/Platypus/main"

workflow {

  input_vcf = Channel.fromPath(params.samples).
      splitCsv(header: true).
      map{row ->
        tuple(row.dataset.toString(), row.patient.toString(), row.sample.toString(), file(row.vcf), file(row.vcf_tbi))}

  input_CNA = Channel.fromPath(params.samples).
    splitCsv(header: true).
    map{row ->
     tuple(row.dataset.toString(), row.patient.toString(), row.sample.toString(), file(row.cna_calling))}

  normal_bam = Channel.fromPath(params.samples).
    splitCsv(header: true).
    map{row ->
     tuple(row.dataset.toString(), row.patient.toString(), row.sample.toString(), file(row.normal_bam), file(row.normal_bai))}

  tumor_bam = Channel.fromPath(params.samples).
    splitCsv(header: true).
    map{row ->
     tuple(row.dataset.toString(), row.patient.toString(), row.sample.toString(), file(row.tumour_bam), file(row.tumour_bai))}  


  //VARIANT_ANNOTATION(input_vcf)
  //PLATYPUS_CALL_VARIANTS(input_vcf.groupTuple(by: [0,1]), normal_bam.groupTuple(by: [0,1,2,3,4]), tumor_bam.groupTuple(by: [0,1]))  
  MPILEUP()
  join_CNAqc = QC(PLATYPUS_CALL_VARIANTS.out.vcf.transpose(by: [2]), input_CNA)
  //SUBCLONAL_DECONVOLUTION(join_CNAqc)

}
