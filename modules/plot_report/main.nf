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
    library(rhdf5)
    library(patchwork)

    format_sampleID = function(sample_ids) {
        stringr::str_replace_all(sample_ids, pattern="^\\[|\\]\$", replacement="") %>% 
        stringr::str_replace_all(pattern="^\\[|\\]\$", replacement="") %>% 
        stringr::str_replace_all(pattern="\\],\\[", replacement=":") %>% 
        strsplit(':') %>% unlist()
    }

    format_patientID = function(patient_ids){
        stringr::str_replace_all(patient_ids, pattern="^\\[|\\]\$", replacement="") %>% 
        strsplit(', ') %>% unlist()
    }

    patientID = format_patientID("$patientID")
    sampleID = format_sampleID("$sampleID")

    cnaqc_data_plot = format_patientID("$cnaqc_data_plot")
    cnaqc_data_qc = format_patientID("$cnaqc_data_qc")

    viber_pdf = format_patientID("$viber_pdf")
    mobster_pdf = format_patientID("$mobster_pdf")
    ctree_mobster_pdf = format_patientID("$ctree_mobster_pdf")


    #tsv_joint = read.table(file = "$joint_table", header = T, sep = "\t")

    pdf(file = "report/$datasetID/final_report.pdf", paper = 'a4')

    # Dataset 
    maftools_oncoplot <- image_read_pdf("$maftools_oncoplot") %>% image_ggplot()
    maftools_summary <- image_read_pdf("$maftools_summary") %>% image_ggplot()
    maftools <- maftools_summary + maftools_oncoplot + patchwork::plot_layout(nrow = 2)
    maftools

    image_read_pdf("$spareSig_plot")
    
    # Sample
    lapply(1:length(patientID), function(p){

      samples_pID = format_patientID(sampleID[[p]])

      lapply(1:length(samples_pID), function(s){
        
        image_read_pdf( format_list(cnaqc_data_plot[[p]])[[s]] ) %>% image_ggplot()
        image_read_pdf( format_list(cnaqc_qc_plot[[p]])[[s]] ) %>% image_ggplot()

        image_read_pdf( format_list(viber_pdf[[p]])[[s]] ) %>% image_ggplot()
        image_read_pdf( format_list(mobster_pdf[[p]])[[s]] ) %>% image_ggplot()
        image_read_pdf( format_list(ctree_mobster_pdf[[p]])[[s]] ) %>% image_ggplot()
        
        #pyclone_fit = read.table(header = T)
        #pyclone_h5 = 
        #plot_pyclone = plot_summary_pyclone(x = tsv_joint,
        #              y = pyclone_fit,
        #              h5_file = pyclone_h5,
        #              d1 = p)
      
      })
    })

    dev.off()

    """


}