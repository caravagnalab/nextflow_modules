#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { VARIANT_ANNOTATION } from "${baseDir}/subworkflows/variant_annotation/main"
include { FORMATTER } from "${baseDir}/subworkflows/formatter/main"
include { LIFTER } from "${baseDir}/subworkflows/lifter/main"
include { DRIVER_ANNOTATION } from "${baseDir}/subworkflows/driver_annotation/main"
include { QC } from "${baseDir}/subworkflows/QC/main"
include { SUBCLONAL_DECONVOLUTION } from "${baseDir}/subworkflows/subclonal_deconvolution/main"


workflow {

  input_vcf = Channel.fromPath(params.samples).
      splitCsv(header: true).
      map{row ->
        tuple(row.dataset.toString(), row.patient.toString(), row.sample.toString(), file(row.vcf), file(row.vcf_tbi))}

  input_CNA = Channel.fromPath(params.samples).
    splitCsv(header: true).
    map{row ->
     tuple(row.dataset.toString(), row.patient.toString(), row.sample.toString(), file(row.cna_calling))}

  // normal_bam = Channel.fromPath(params.samples).
  //   splitCsv(header: true).
  //   map{row ->
  //    tuple(row.dataset.toString(), row.patient.toString(), row.sample.toString(), file(row.normal_bam), file(row.normal_bai))}

  tumor_bam = Channel.fromPath(params.samples).
    splitCsv(header: true).
    map{row ->
     tuple(row.dataset.toString(), row.patient.toString(), row.sample.toString(), file(row.tumour_bam), file(row.tumour_bai))}  


  VARIANT_ANNOTATION(input_vcf)
  FORMATTER(VARIANT_ANNOTATION.out.vep_out)
  FORMATTER(input_CNA)

  // if multisample
  LIFTER(FORMATTER.out, tumor_bam)
  DRIVER_ANNOTATION(LIFTER.out)

  // if singlesample
  DRIVER_ANNOTATION(FORMATTER.out.vcf)

  // both
  join_CNAqc = QC(DRIVER_ANNOTATION.out, FORMATTER.out.cna)

  //SUBCLONAL_DECONVOLUTION(join_CNAqc)

}
