#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { VARIANT_ANNOTATION } from "${baseDir}/subworkflows/variant_annotation/main"
include { SUBCLONAL_DECONVOLUTION } from "${baseDir}/subworkflows/subclonal_deconvolution/main"
include { CNAQC } from "${baseDir}/subworkflows/CNAqc/main"

workflow {

  input_vcf = Channel.fromPath(params.samples).
      splitCsv(header: true).
      map{row ->
        tuple(row.dataset.toString(), row.patient.toString(), row.sample.toString(), file(row.vcf))}

  input_CNA = Channel.fromPath(params.samples).
    splitCsv(header: true).
    map{row ->
     tuple(row.dataset.toString(), row.patient.toString(), row.sample.toString(), file(row.cna_calling))}

  annotated_vcf = VARIANT_ANNOTATION(input_vcf)
  join_CNAqc = CNAQC(annotated_vcf, input_CNA)
  //SUBCLONAL_DECONVOLUTION(join_CNAqc)

}
