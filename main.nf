#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { SUBCLONAL_DECONVOLUTION } from "${baseDir}/subworkflows/subclonal_deconvolution/main"
include { VARIANT_ANNOTATION } from "${baseDir}/subworkflows/variant_annotation/main"


workflow {

  input_vcf = Channel.fromPath(params.samples).
      splitCsv(header: true).
      map{row ->
        tuple(row.dataset.toString(), row.patient.toString(), row.sample.toString(), file(row.vcf))}
  

  input_sequenza = Channel.fromPath(params.samples).
    splitCsv(header: true).
    map{row ->
     tuple(row.dataset.toString(), row.patient.toString(), row.sample.toString(), row.sex.toString(), file(row.seqz), file(row.vcf))}

annotation_output = ANNOVAR_ANNOTATE(input_vcf)
vcf2maf_out = VCF2MAF(vep_out)
maf_out = MAFTOOLS(vcf2maf_out.groupTuple(by: 0))
annovar2maf_out = ANNOVAR2MAF(annovar_out)
SUBCLONAL_DECONVOLUTION(joint_table)
}
