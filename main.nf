#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SEQUENZA_EXTRACT } from "/u/cdslab/ncalonaci/sequenza_tests/sequenza/main.nf"

Channel
    .fromPath( params.seqz )
    .ifEmpty { error "Cannot find any reads matching: ${params.seqz}" }
    .set(seqz) 

workflow{
  SEQUENZA_EXTRACT(seqz)
}

