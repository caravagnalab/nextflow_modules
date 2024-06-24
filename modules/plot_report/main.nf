process PLOT_REPORT_SINGLE_SAMPLE {
  publishDir params.publish_dir, mode: 'copy'
  
  input:
    
    tuple val(dID_maftool), val(pID_maftool), val(sID_maftool), path(maftools_oncoplot) 
    tuple val(dID_maftool), val(pID_maftool), val(sID_maftool), path(maftools_summary)
    tuple val(dID_cnaqcData), val(pID_cnaqcData), val(sID_cnaqcData), path(cnaqc_data_plot,  stageAs: 'cnaqc_data*.pdf')
    tuple val(dID_cnaqcQC), val(pID_cnaqcQC), val(sID_cnaqcQC), path(cnaqc_qc_plot,  stageAs: 'cnaqc_qc*.pdf')
    tuple val(dID_sig), val(pID_sig), val(sID_sig), path(spareSig_plot)
    tuple val(dID_viber), val(pID_viber), val(sID_viber), path(viber_pdf, stageAs: 'viber*.pdf')
    //tuple val(datasetID),   val(patientID), val(sampleID), path(ctree_viber_pdf)
    tuple val(dID_pycloneFit), val(pID_pycloneFit), val(sID_pycloneFit), path(pyclone_fits, stageAs: 'pyclone*.h5')
    tuple val(dID_pycloneBest), val(pID_pycloneBest), val(sID_pycloneBest), path(pyclone_best, stageAs: 'pyclone_best*.txt')
    //tuple val(datasetID),   val(patientID), val(sampleID), path(ctree_pyclone_pdf)
    tuple val(dID_mobster), val(pID_mobster), val(sID_mobster), path(mobster_pdf, stageAs: 'mobster*.pdf')
    tuple val(dID_ctreeMob), val(pID_ctreeMob), val(sID_ctreeMob), path(ctree_mobster_pdf, stageAs: 'ctree_mob*.pdf')
    tuple val(dID_table), val(pID_table), val(sID_table), path(pyclone_table, stageAs: 'pyclone_table*.tsv')
  
  
  output:

    path("report/$dID_sig/final_report.pdf"), emit: pdf

  script:

    """
    #!/usr/bin/env Rscript
    source(paste0("$moduleDir", '/pyclone_plot.R'))

    library(tidyverse)
    library(ggplot2)
    library(magick)
    library(rhdf5)
    library(patchwork)

    format_ID = function(ids) {
      stringr::str_replace_all(ids, pattern="^\\\\[|\\\\]\$", replacement="") %>% 
      strsplit(', ') %>% unlist()
    }

    split_path = function(path){
        strsplit(path, ' ') %>% unlist()
    }

    cnaqc_data_plot = split_path("$cnaqc_data_plot")
    names(cnaqc_data_plot) = format_ID("$sID_cnaqcData")

    cnaqc_qc_plot = split_path("$cnaqc_qc_plot")
    names(cnaqc_qc_plot) = format_ID("$sID_cnaqcQC")

    table_pyclone = split_path("$pyclone_table")
    names(table_pyclone) = format_ID("$sID_table")

    pyclone_fits = split_path("$pyclone_fits")
    names(pyclone_fits) = format_ID("$sID_pycloneFit")

    pyclone_best = split_path("$pyclone_best")
    names(pyclone_best) = format_ID("$sID_pycloneBest")

    viber_pdf = split_path("$viber_pdf")
    names(viber_pdf) = format_ID("$sID_viber")

    mobster_pdf = split_path("$mobster_pdf")
    names(mobster_pdf) = format_ID("$sID_mobster")

    ctree_mobster_pdf = split_path("$ctree_mobster_pdf")
    names(ctree_mobster_pdf) = format_ID("$sID_ctreeMob")
    
    res_dir = paste0("report/", "$dID_sig", "/")
    dir.create(res_dir, recursive = T)
    pdf(file = paste0(res_dir, "final_report.pdf"), paper = 'a4')

    # Dataset 
    maftools_oncoplot = image_read_pdf("$maftools_oncoplot") %>% image_ggplot()
    maftools_summary = image_read_pdf("$maftools_summary") %>% image_ggplot()
    maftools = maftools_summary + maftools_oncoplot + patchwork::plot_layout(nrow = 2)
    print(maftools)

    signatures = image_read_pdf("$spareSig_plot") %>% image_ggplot()
    print(signatures)
    
    sampleID = format_ID("$sID_cnaqcData")
    for (i in sampleID){
      print(image_read_pdf(cnaqc_data_plot[[i]]) %>% image_ggplot())
      print(image_read_pdf(cnaqc_qc_plot[[i]]) %>% image_ggplot())
      print(image_read_pdf(viber_pdf[[i]]) %>% image_ggplot())
      print(image_read_pdf(mobster_pdf[[i]]) %>% image_ggplot())
      print(image_read_pdf(ctree_mobster_pdf[[i]]) %>% image_ggplot())

      # PLOT PYCLONE
      table =  read.table(file = table_pyclone[[i]], header = T, sep = "\\t")
      best = read.table(pyclone_best[[i]], header = T)
      fit = pyclone_fits[[i]]
      plot_pyclone = plot_summary_pyclone(
                      x = table,
                      y = best,
                      h5_file = fit,
                      d1 = i)
      print(plot_pyclone)
    }

    # Sample
    dev.off()

    """


}
