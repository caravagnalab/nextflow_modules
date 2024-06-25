#!/usr/bin/env nextflow
nextflow.enable.dsl=2


include { VARIANT_ANNOTATION } from "${baseDir}/subworkflows/variant_annotation/main"
include { FORMATTER as FORMATTER_CNA } from "${baseDir}/subworkflows/formatter/main"
include { FORMATTER as FORMATTER_VCF} from "${baseDir}/subworkflows/formatter/main"
include { LIFTER } from "${baseDir}/subworkflows/lifter/main"
include { DRIVER_ANNOTATION } from "${baseDir}/subworkflows/annotate_driver/main"
include { FORMATTER as FORMATTER_RDS} from "${baseDir}/subworkflows/formatter/main"
include { QC } from "${baseDir}/subworkflows/QC/main"
include { SUBCLONAL_DECONVOLUTION } from "${baseDir}/subworkflows/subclonal_deconvolution/main"
include { MUTATIONAL_SIGNATURES } from "${baseDir}/subworkflows/mutational_signatures/main"
include { PLOT_REPORT_SINGLE_SAMPLE } from "${baseDir}/modules/plot_report/main"
include { PLOT_REPORT_MULTI_SAMPLE } from "${baseDir}/modules/plot_report/plot_report_multi"


workflow {  

  input_vcf = Channel.fromPath(params.samples).
      splitCsv(header: true).
      map{row ->
        tuple(row.dataset.toString(), row.patient.toString(), row.sample.toString(), file(row.vcf), file(row.vcf_tbi))}
 
  cancer_type = Channel.fromPath(params.samples).
      splitCsv(header: true).
      map{row -> row.cancer_type.toString()}
 
  input_cna = Channel.fromPath(params.samples).
    splitCsv(header: true).
    map{row ->
     tuple(row.dataset.toString(), row.patient.toString(), row.sample.toString(), file(row.cnv_res), row.cnv_caller.toString())}

  VARIANT_ANNOTATION(input_vcf)
  FORMATTER_VCF(VARIANT_ANNOTATION.out.vep, "vcf")
  FORMATTER_CNA(input_cna, "cna")
  
  exist_bam_val = false //placeholder
  if (params.mode == 'multisample' && exist_bam_val){  
    tumor_bam = Channel.fromPath(params.samples).
      splitCsv(header: true).
       map{row ->
       tuple(row.dataset.toString(), row.patient.toString(), row.sample.toString(), file(row.tumour_bam), file(row.tumour_bai))}
    
    LIFTER(FORMATTER_VCF.out, tumor_bam)
    annotation = DRIVER_ANNOTATION(LIFTER.out, cancer_type)

  } else {
   annotation = DRIVER_ANNOTATION(FORMATTER_VCF.out, cancer_type)
  }
  
  QC(FORMATTER_CNA.out, annotation)
  SUBCLONAL_DECONVOLUTION(QC.out.rds_join)
  MUTATIONAL_SIGNATURES(QC.out.rds_join)

  if (params.mode == 'singlesample'){
     PLOT_REPORT_SINGLE_SAMPLE(
             VARIANT_ANNOTATION.out.maf_oncoplot, 
             VARIANT_ANNOTATION.out.maf_summary_plot,
             QC.out.plot_cnaqc_data.groupTuple(by: 0), 
             QC.out.plot_cnaqc_qc.groupTuple(by: 0),
             MUTATIONAL_SIGNATURES.out.plot_pdf,
             SUBCLONAL_DECONVOLUTION.out.viber_pdf.groupTuple(by: 0),
             SUBCLONAL_DECONVOLUTION.out.pyclone_fits.groupTuple(by: 0),
             SUBCLONAL_DECONVOLUTION.out.pyclone_best.groupTuple(by: 0),
             SUBCLONAL_DECONVOLUTION.out.mobster_pdf.groupTuple(by: 0),
             SUBCLONAL_DECONVOLUTION.out.ctree_mobster_pdf.groupTuple(by: 0),
             SUBCLONAL_DECONVOLUTION.out.pyclone_table.groupTuple(by: 0)
             )
  } else if  (params.mode == 'multisample') {
      PLOT_REPORT_MULTI_SAMPLE(
            VARIANT_ANNOTATION.out.maf_oncoplot, 
            VARIANT_ANNOTATION.out.maf_summary_plot,
            QC.out.plot_cnaqc_data.groupTuple(by: 0), 
            QC.out.plot_cnaqc_qc.groupTuple(by: 0),
            MUTATIONAL_SIGNATURES.out.plot_pdf,
            SUBCLONAL_DECONVOLUTION.out.viber_pdf.groupTuple(by: 0),
            SUBCLONAL_DECONVOLUTION.out.ctree_viber_pdf.groupTuple(by: 0),
            SUBCLONAL_DECONVOLUTION.out.pyclone_fits.groupTuple(by: 0),
            SUBCLONAL_DECONVOLUTION.out.pyclone_best.groupTuple(by: 0),
            //SUBCLONAL_DECONVOLUTION.out.ctree_pyclone_pdf.groupTuple(by: [0]),
            SUBCLONAL_DECONVOLUTION.out.pyclone_table.groupTuple(by: 0)
            )
  }
}

