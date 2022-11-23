#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { SEQUENZA_EXTRACT } from "${baseDir}/sequenza/main"
// include { PARSE_INPUT } from "${baseDir}/sequenza/main"

workflow {
  // WORKING
  // Channel.fromPath(params.samples) \
  //   | splitCsv(header: true) \
  //   | map{row ->
  //       tuple(row.patient.toString(), row.sample.toString(), file(row.file))} \
  //   | SEQUENZA_EXTRACT
  
  input_rows = Channel.fromPath(params.samples).
    splitCsv(header: true).
    map{row ->
      tuple(row.patient.toString(), row.sample.toString(), file(row.file))}  
  SEQUENZA_EXTRACT(input_rows)
}
