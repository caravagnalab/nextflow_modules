#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Channel
//     .fromPath( params.seqz )
//     .ifEmpty { error "Cannot find any reads matching: ${params.seqz}" } // if empty, complains
//     .set {seqzFile}

include { SEQUENZA_EXTRACT } from "${baseDir}/sequenza/main" addParams(OUTPUT: params.out_dir)

workflow {
  SEQUENZA_EXTRACT(params.seqz)
}
