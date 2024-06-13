process PLOT_REPORT_MULTI_SAMPLE {
  publishDir params.publish_dir, mode: 'copy'
  
  input:
    
    tuple val(datasetID),   val(patientID), val(sampleID), path(maftools_oncoplot) 
    tuple val(datasetID),   val(patientID), val(sampleID), path(maftools_summary)
    tuple val(datasetID),   val(patientID), val(sampleID), path(cnaqc_data_plot,  stageAs: 'cnaqc_data*.pdf')
    tuple val(datasetID),   val(patientID), val(sampleID), path(cnaqc_qc_plot,  stageAs: 'cnaqc_qc*.pdf')
    tuple val(datasetID),   val(patientID), val(sampleID), path(spareSig_plot)
    tuple val(datasetID),   val(patientID), val(sampleID), path(viber_pdf)
    tuple val(datasetID),   val(patientID), val(sampleID), path(ctree_viber_pdf)
    tuple val(datasetID),   val(patientID), val(sampleID), path(pyclone_fits)
    tuple val(datasetID),   val(patientID), val(sampleID), path(pyclone_best)
    tuple val(datasetID),   val(patientID), val(sampleID), path(ctree_pyclone_pdf)
    //tuple val(datasetID),   val(patientID), val(sampleID), path(mobster_pdf)

  
  output:

    tuple val(datasetID), val(patientID), val (sampleID), path("report/$datasetID/final_report.pdf"), emit: pdf

  script:

    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    library(ggplot2)
    library(magick)

    patientID = llll"$patientID"
    sampleID = llll"$sampleID"

    cnaqc_data_plot = llll"$cnaqc_data_plot"
    cnaqc_data_qc = llll"$cnaqc_data_qc"



    open(pdf)

    # Dataset 
    maftools_oncoplot <- image_read_pdf("$maftools_oncoplot") %>% image_ggplot()
    maftools_summary <- image_read_pdf("$maftools_summary") %>% image_ggplot()
    maftools <- maftools_oncoplot | maftools_summary
    maftools

    image_read_pdf("$spareSig_plot")
    
    # Patient

    # Sample
    lapply(1:length(patientID), function(p){

      lapply(1:length(sampleID[[p]]), function(s){
        
        image_read_pdf(cnaqc_data_plot[[p]][[s]]) 
        image_read_pdf(cnaqc_qc_plot[[p]][[s]]) 


      
      })
    })

  


    close(pdf)






    """


}