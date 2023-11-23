#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//include { SEQUENZA_EXTRACT } from "${baseDir}/nextflow_modules/modules/sequenza/main"
include { VEP_ANNOTATE } from "${baseDir}/nextflow_modules/modules/VEP/main"
include { VCF2MAF } from "${baseDir}/nextflow_modules/modules/vcf2maf/main"
include { MAFTOOLS } from "${baseDir}/nextflow_modules/modules/maftools/main_vep"
//include { ANNOVAR_ANNOTATE } from "${baseDir}/nextflow_modules/modules/annovar/main"
//include { ANNOVAR2MAF } from "${baseDir}/nextflow_modules/modules/annovar2maf/main2" 
//include { SNPEFF_ANNOTATE } from "${baseDir}/nextflow_modules/modules/snpeff/main" 
//include { BCFTOOLS_SPLIT_VEP } from "${baseDir}/nextflow_modules/modules/bcftools/main"
//PLATYPUS_CALL_VARIANTS(input_multisample.groupTuple(by: [0,3,4]))include { PLATYPUS_CALL_VARIANTS } from "${baseDir}/modules/Platypus/main"
//include { JOIN_TABLES } from "${baseDir}/modules/join_tables/main"
//include { SEQUENZA_CNAqc } from "${baseDir}/modules/Sequenza_CNAqc/main"

workflow {

  input_vcf = Channel.fromPath(params.samples).
      splitCsv(header: true).
      map{row ->
        tuple(row.patient.toString(), row.sample.toString(), file(row.vcf))}
  
}

  //input_bam = Channel.fromPath(params.samples).
  //splitCsv(header: true).
  //  map{row ->
  //    tuple(row.patient.toString(), row.timepoint.toString(), row.sample.toString(), file(row.tumor_bam), file(row.nor//mal_bam))}

  //input_sequenza = Channel.fromPath(params.samples).
  //  splitCsv(header: true).
  //  map{row ->
  //    tuple(row.patient.toString(), row.timepoint.toString(), row.sample.toString(), row.sex.toString(), file(row.seqz//_file), file(row.snv_vcf_f//ile))}

  //input_multisample = Channel.fromPath(params.samples).
  //    splitCsv(header: true).
  //    map{row ->
  //      tuple(row.patient.toString(), file(row.tumor_bam), file(row.tumor_bai), file(row.normal_bam), file(row.normal_//bai), file(row.snv_vcf_fil//e), file(row.snv_tbi_file), file(row.indel_vcf_file), file(row.indel_tbi_file))}

vep_output = VEP_ANNOTATE(input_vcf)
  //BCFTOOLS_SPLIT_VEP(vep_output)
  //PLATYPUS_CALL_VARIANTS(input_multisample.groupTuple(by: [0,3,4]))
  //JOINT_TABLE()
  //SEQUENZA_CNAqc(input_sequenza)
vcf2maf_output = VCF2MAF(vep_output)
maf_output = MAFTOOLS(vcf2maf_output.groupTuple()) 
//annovar_output = ANNOVAR_ANNOTATE(input_vcf)
//annovar2maf_output = ANNOVAR2MAF(annovar_output)
//snpeff_output = SNPEFF_ANNOTATE(input_vcf)
 
  
}
                     
