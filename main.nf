#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// include { MOBSTERh } from "${baseDir}/modules/mobsterh/main"
// include { PYCLONEVI } from "${baseDir}/modules/pyclonevi/main"
// include { VIBER } from "${baseDir}/modules/viber/main"
// include { CTREE } from "${baseDir}/modules/ctree/main"
// include { VARTRIX } from "${baseDir}/modules/vartrix/main"
include { SUBCLONAL_DECONVOLUTION } from "${baseDir}/subworkflows/subclonal_deconvolution/main"

workflow {

  input_subclonal = Channel.fromPath(params.samples).
       splitCsv(header: true).
       map{row -> tuple(row.patient.toString(), file(row.joint_table))}

  SUBCLONAL_DECONVOLUTION(input_subclonal)

  // mobster = MOBSTERh(input_subclonal)
  // CTREE(mobster)
  // vartrix = VARTRIX(input_vartrix)
  // pyclonevi = PYCLONEVI(input_pyclonevi)

  }
