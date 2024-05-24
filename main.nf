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
include { FORMATTER as FORMATTER_RDS} from "${baseDir}/subworkflows/formatter/main"
include { SUBCLONAL_DECONVOLUTION } from "${baseDir}/subworkflows/subclonal_deconvolution/main"


workflow {

  //input_vcf = Channel.fromPath(params.samples).
  //    splitCsv(header: true).
  //    map{row ->
  //      tuple(row.dataset.toString(), row.patient.toString(), row.sample.toString(), file(row.vcf))}
  

  //input_sequenza = Channel.fromPath(params.samples).
  //  splitCsv(header: true).
  //  map{row ->
  //   tuple(row.patient.toString(), row.sample.toString(), row.sex.toString(), file(row.seqz), file(row.vcf))

  input_joint_table = Channel.fromPath(params.samples).
            splitCsv(header: true).
            map{row ->
              tuple(row.dataset.toString(), row.patient.toString(), row.sample.toString(),file(row.joint_table))}.groupTuple(by: [0,1,3])

//vep_output = VEP_ANNOTATE(input_vcf) 
//annovar_output = ANNOVAR_ANNOTATE(input_vcf)
//vcf2maf_output = VCF2MAF(vep_output)
//maf_output = MAFTOOLS(vcf2maf_output.groupTuple(by: 0))
//BCFTOOLS_SPLIT_VEP(vep_output)
//PLATYPUS_CALL_VARIANTS(input_multisample.groupTuple(by: [0,3,4]))
//JOINT_TABLE()
//SEQUENZA_CNAqc(input_sequenza)
FORMATTER_RDS(input_joint_table, "rds")
//SUBCLONAL_DECONVOLUTION(FORMATTER_RDS.out)
}
