#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//include { SEQUENZA_EXTRACT } from "${baseDir}/modules/sequenza/main"
// include { VEP_ANNOTATE } from "${baseDir}/modules/VEP/main"
include { MOBSTERh } from "${baseDir}/modules/mobsterh/main"
//include { VCF2MAF } from "${baseDir}/modules/vcf2maf/main"
//include { BCFTOOLS_SPLIT_VEP } from "${baseDir}/modules/bcftools/main"
//include { SEQUENZA_CNAqc } from "${baseDir}/modules/Sequenza_CNAqc/main"
//include { JOIN_TABLES } from "${baseDir}/modules/join_tables/main"

workflow {

  // input_vcf = Channel.fromPath(params.samples).
  //     splitCsv(header: true).
  //     map{row ->
  //       tuple(row.patient.toString(), row.sample.toString(), file(row.vcf))} 

  //input_sequenza = Channel.fromPath(params.samples).
  //  splitCsv(header: true).
  //  map{row ->
  //    tuple(row.patient.toString(), row.sample.toString(), row.sex.toString(), file(row.seqz), file(row.vcf))}

  //input_multisample = Channel.fromPath(params.samples).
  //    splitCsv(header: true).
  //    map{row ->
  //      tuple(row.patient.toString(), file(row.tumour_bam), file(row.tumour_bai), file(row.normal_bam), file(row.normal_bai), file(row.vcf), file(row.vcf_tbi))} 

  input_mobsterh = Channel.fromPath(params.samples).
        splitCsv(header: true).
        map{row ->
          tuple(row.patient.toString(), row.timepoint.toString(), row.sample.toString(), file(row.joint_table))} 

  mobster = MOBSTERh(input_mobsterh)

  // vep_output = VEP_ANNOTATE(input_vcf)
  //VCF2MAF(input_vcf)
  //BCFTOOLS_SPLIT_VEP(vep_output)
  //PLATYPUS_CALL_VARIANTS(input_multisample.groupTuple(by: [0,3,4]))
  //JOINT_TABLE()
  //SEQUENZA_CNAqc(input_sequenza)
  }
