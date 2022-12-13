#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { SEQUENZA_EXTRACT } from "${baseDir}/sequenza/main"
include { VEP_ANNOTATE } from "${baseDir}/VEP/main"
// include { PARSE_INPUT } from "${baseDir}/sequenza/main"

workflow {
  // WORKING
  // Channel.fromPath(params.samples) \
  //   | splitCsv(header: true) \
  //   | map{row ->
  //       tuple(row.patient.toString(), row.sample.toString(), file(row.file))} \
  //   | SEQUENZA_EXTRACT
  
  input_sequenza = Channel.fromPath(params.samples).
    splitCsv(header: true).
    map{row ->
      tuple(row.patient.toString(), row.sample.toString(), file(row.seqz_file))}  
  input_vcf = Channel.fromPath(params.samples).
      splitCsv(header: true).
      map{row ->
        tuple(row.patient.toString(), row.sample.toString(), file(row.vcf_file))} 
  SEQUENZA_EXTRACT(input_sequenza)
  VEP_ANNOTATE(input_vcf)
}
