#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//include { VEP_ANNOTATE } from "${baseDir}/modules/VEP/main"
//include { VCF2MAF } from "${baseDir}/modules/vcf2maf/main"
//include { MAFTOOLS } from "${baseDir}/modules/maftools/main"
//include { ANNOVAR_ANNOTATE } from "${baseDir}/modules/annovar/main"
//include { SEQUENZA_CNAqc } from "${baseDir}/modules/Sequenza_CNAqc/main"
//include { BCFTOOLS_SPLIT_VEP } from "${baseDir}/modules/bcftools/main"
//include { JOIN_TABLES } from "${baseDir}/modules/join_tables/main"
//include { SEQUENZA_CNAqc } from "${baseDir}/modules/Sequenza_CNAqc/main"
include { VARIANT_ANNOTATION } from "${baseDir}/subworkflows/variant_annotation/main"
include { FORMATTER as FORMATTER_CNA } from "${baseDir}/subworkflows/formatter/main"
include { FORMATTER as FORMATTER_VCF} from "${baseDir}/subworkflows/formatter/main"
include { FORMATTER as FORMATTER_RDS} from "${baseDir}/subworkflows/formatter/main"
include { QC } from "${baseDir}/subworkflows/QC/main"
include { SUBCLONAL_DECONVOLUTION } from "${baseDir}/subworkflows/subclonal_deconvolution/main"

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
  FORMATTER_VCF(VARIANT_ANNOTATION.out.vep, "vcf")//VARIANT_ANNOTATION.out.vep
  FORMATTER_CNA(input_cna, "cna")

  // // if multisample 
  // LIFTER(FORMATTER_VCF.out, tumor_bam) // work
  // DRIVER_ANNOTATION(LIFTER.out)

  // // if singlesample 
  // DRIVER_ANNOTATION(FORMATTER.out.vcf)

  join_CNAqc = QC(FORMATTER_CNA.out, FORMATTER_VCF.out) // work

  SUBCLONAL_DECONVOLUTION(join_CNAqc)

}
