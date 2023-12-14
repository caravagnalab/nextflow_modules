#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { MOBSTERh } from "${baseDir}/modules/mobsterh/main"
// include { PYCLONEVI } from "${baseDir}/modules/pyclonevi/main"
// include { VIBER } from "${baseDir}/modules/viber/main"
include { CTREE } from "${baseDir}/modules/ctree/main"
// include { VARTRIX } from "${baseDir}/modules/vartrix/main"

workflow {

  input_mobsterh = Channel.fromPath(params.samples).
       splitCsv(header: true).
       map{row ->
         tuple(row.patient.toString(), row.timepoint.toString(), row.sample.toString(), file(row.joint_table))}

  //input_vartrix = Channel.fromPath(params.samples).
  //    splitCsv(header: true).
  //    map{row ->
  //      tuple(row.patient.toString(), row.sample.toString(), file(row.bam), file(row.bai), file(row.vcf), file(row.barcode))}

  // input_pyclonevi = Channel.fromPath(params.samples).
  //     splitCsv(header: true).
  //     map{row ->
  //       tuple(row.patient.toString(), file(row.joint_table))}

  mobster = MOBSTERh(input_mobsterh)
  CTREE(mobster)
  // vartrix = VARTRIX(input_vartrix)
  // pyclonevi = PYCLONEVI(input_pyclonevi)
  // viber = VIBER(input_pyclonevi)

  }
