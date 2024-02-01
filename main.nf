#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { VARIANT_ANNOTATION } from "${baseDir}/subworkflows/variant_annotation/main"
include { SUBCLONAL_DECONVOLUTION } from "${baseDir}/subworkflows/subclonal_deconvolution/main"

workflow {

  input_vcf = Channel.fromPath(params.samples).
      splitCsv(header: true).
      map{row ->
        tuple(row.dataset.toString(), row.patient.toString(), row.sample.toString(), file(row.vcf))}
  

  input_sequenza = Channel.fromPath(params.samples).
    splitCsv(header: true).
    map{row ->
     tuple(row.dataset.toString(), row.patient.toString(), row.sample.toString(), row.sex.toString(), file(row.seqz), file(row.vcf))}

annotated_vcf = VARIANT_ANNOTATION(input_vcf)
//SUBCLONAL_DECONVOLUTION(joint_table)
}
