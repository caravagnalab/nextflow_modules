process PLOT_REPORT_SINGLE_SAMPLE {
  publishDir params.publish_dir, mode: 'copy'
  
  input:
    
    tuple val(datasetID),   val(patientID), val(sampleID), path(maftools_oncoplot) 
    tuple val(datasetID),   val(patientID), val(sampleID), path(maftools_summary)
    tuple val(datasetID),   val(patientID), val(sampleID), path(cnaqc_data_plot,  stageAs: 'cnaqc_data*.pdf')
    tuple val(datasetID),   val(patientID), val(sampleID), path(cnaqc_qc_plot,  stageAs: 'cnaqc_qc*.pdf')
    tuple val(datasetID),   val(patientID), val(sampleID), path(spareSig_plot)
    tuple val(datasetID),   val(patientID), val(sampleID), path(viber_pdf)
    //tuple val(datasetID),   val(patientID), val(sampleID), path(ctree_viber_pdf)
    tuple val(datasetID),   val(patientID), val(sampleID), path(pyclone_fits)
    tuple val(datasetID),   val(patientID), val(sampleID), path(pyclone_best)
    //tuple val(datasetID),   val(patientID), val(sampleID), path(ctree_pyclone_pdf)
    tuple val(datasetID),   val(patientID), val(sampleID), path(mobster_pdf)
    tuple val(datasetID),   val(patientID), val(sampleID), path(ctree_mobster_pdf)
  
  output:

    tuple val(datasetID), val(patientID), val (sampleID), path("report/$datasetID/final_report.pdf"), emit: pdf

  script:

    """
    #!/usr/bin/env Rscript
    source(paste0("$moduleDir", '/pyclone_plot.R'))

    library(tidyverse)
    library(ggplot2)
    library(magick)

    format_list = function(nf_list) {
      stringr::str_replace_all(nf_list, pattern="^\\[|\\]\$", replacement="") %>% 
        stringr::str_replace_all(pattern="\\], \\[", replacement="\\]:\\[") %>% 
        strsplit(':') %>% unlist()
    }

    patientID = format_list("$patientID")
    sampleID = format_list("$sampleID")

    cnaqc_data_plot = format_list("$cnaqc_data_plot")
    cnaqc_data_qc = format_list("$cnaqc_data_qc")

    viber_pdf = format_list("$viber_pdf")
    mobster_pdf = format_list("$mobster_pdf")
    ctree_mobster_pdf = format_list("$ctree_mobster_pdf")


    # tsv_joint <- read.table(file = "~/single_sample/CNAqc2tsv/TEST/MSeq_Set06/joint_table.tsv",header = T, sep = "\t")
    # pyclone_fit <- read.table(file = "~/MSeq_Set06/best_fit.txt", header = T)
    # plot_pyclone = plot_summary_pyclone(x = tsv_joint,
    #                       y = pyclone_fit,
    #                       h5_file = "/path/to/all_fits.h5",
    #                       d1 = "Set6_42",
    #                       d2 = "Set6_44")

    pdf(file = "report/$datasetID/final_report.pdf", paper = 'a4')

    # Dataset 
    maftools_oncoplot <- image_read_pdf("$maftools_oncoplot") %>% image_ggplot()
    maftools_summary <- image_read_pdf("$maftools_summary") %>% image_ggplot()
    maftools <- maftools_oncoplot | maftools_summary
    maftools

    image_read_pdf("$spareSig_plot")
    
    # Patient

    # Sample
    lapply(1:length(patientID), function(p){

      samples_pID = format_list(sampleID[[p]])

      lapply(1:length(samples_pID), function(s){
        
        image_read_pdf( format_list(cnaqc_data_plot[[p]])[[s]] ) %>% image_ggplot()
        image_read_pdf( format_list(cnaqc_qc_plot[[p]])[[s]] ) %>% image_ggplot()

        image_read_pdf( format_list(viber_pdf[[p]])[[s]] ) %>% image_ggplot()
        image_read_pdf( format_list(mobster_pdf[[p]])[[s]] ) %>% image_ggplot()
        image_read_pdf( format_list(ctree_mobster_pdf[[p]])[[s]] ) %>% image_ggplot()
        # plot_pyclone
      
      })
    })

    dev.off()






    """


}