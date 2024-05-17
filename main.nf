#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { VARIANT_ANNOTATION } from "${baseDir}/subworkflows/variant_annotation/main"
include { FORMATTER as FORMATTER_CNA } from "${baseDir}/subworkflows/formatter/main"
include { FORMATTER as FORMATTER_VCF} from "${baseDir}/subworkflows/formatter/main"
// include { LIFTER } from "${baseDir}/subworkflows/lifter/main"
// include { DRIVER_ANNOTATION } from "${baseDir}/subworkflows/driver_annotation/main"
include { QC } from "${baseDir}/subworkflows/QC/main"
// include { SUBCLONAL_DECONVOLUTION } from "${baseDir}/subworkflows/subclonal_deconvolution/main"


workflow {

  input_vcf = Channel.fromPath(params.samples).
      splitCsv(header: true).
      map{row ->
        tuple(row.dataset.toString(), row.patient.toString(), row.sample.toString(), file(row.vcf), file(row.vcf_tbi))}
 
  input_cna = Channel.fromPath(params.samples).
    splitCsv(header: true).
    map{row ->
     tuple(row.dataset.toString(), row.patient.toString(), row.sample.toString(), file(row.cnv_res), row.cnv_caller.toString())}

  // normal_bam = Channel.fromPath(params.samples).
  //   splitCsv(header: true).
  //   map{row ->
  //    tuple(row.dataset.toString(), row.patient.toString(), row.sample.toString(), file(row.normal_bam), file(row.normal_bai))}

  // tumor_bam = Channel.fromPath(params.samples).
  //   splitCsv(header: true).
  //   map{row ->
  //    tuple(row.dataset.toString(), row.patient.toString(), row.sample.toString(), file(row.tumour_bam), file(row.tumour_bai))}  


  VARIANT_ANNOTATION(input_vcf) //work
  FORMATTER_VCF(VARIANT_ANNOTATION.out.vep_out , "vcf")
  FORMATTER_CNA(input_cna, "cna")

// // if multisample 
//   LIFTER(FORMATTER_VCF.out, tumor_bam) // work
//   DRIVER_ANNOTATION(LIFTER.out)

// // if singlesample 
//   DRIVER_ANNOTATION(FORMATTER.out.vcf)

  join_CNAqc = QC(FORMATTER_VCF.out, FORMATTER_CNA.out) // work

  //  SUBCLONAL_DECONVOLUTION(join_CNAqc)

}
