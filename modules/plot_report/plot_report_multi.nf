process PLOT_REPORT_MULTI_SAMPLE {
  publishDir params.publish_dir, mode: 'copy'
  
  input:
    
    tuple val(dID_maftool), val(pID_maftool), val(sID_maftool), path(maftools_oncoplot) 
    tuple val(dID_maftool), val(pID_maftool), val(sID_maftool), path(maftools_summary)
    tuple val(dID_cnaqcData), val(pID_cnaqcData), val(sID_cnaqcData), path(cnaqc_data_plot,  stageAs: 'cnaqc_data*.pdf')
    tuple val(dID_cnaqcQC), val(pID_cnaqcQC), val(sID_cnaqcQC), path(cnaqc_qc_plot,  stageAs: 'cnaqc_qc*.pdf')
    tuple val(dID_sig), val(pID_sig), val(sID_sig), path(spareSig_plot)
    tuple val(dID_viber), val(pID_viber), val(sID_viber), path(viber_pdf, stageAs: 'viber*.pdf')
    tuple val(dID_ctreeViber),   val(pID_ctreeViber), val(sID_ctreeViber), path(ctree_viber_pdf, stageAs: 'ctree_viber*.pdf')
    tuple val(dID_pycloneFit), val(pID_pycloneFit), val(sID_pycloneFit), path(pyclone_fits, stageAs: 'pyclone*.h5')
    tuple val(dID_pycloneBest), val(pID_pycloneBest), val(sID_pycloneBest), path(pyclone_best, stageAs: 'pyclone_best*.txt')
    //tuple val(dID_ctreePyclone),   val(pID_ctreePyclone), val(sID_ctreePyclone), path(ctree_pyclone_pdf, stageAs: 'ctree_pyclone*.txt')
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

    patient_ID = format_ID("$pID_sig")
    sID = stringr::str_replace_all("$sID_sig", pattern="^\\\\[|\\\\]\$", replacement="") %>% 
              stringr::str_replace_all(pattern="\\\\], \\\\[", replacement=";") %>% 
              strsplit(';')  %>%  unlist()
  
    sample_ID = lapply(sID, FUN = function(p){
      stringr::str_replace_all(p, pattern="^\\\\[|\\\\]\$", replacement="") %>% 
        strsplit(', ') %>% unlist()
    })
    names(sample_ID) = patient_ID
 
    cnaqc_data_plot = split_path("$cnaqc_data_plot")
    names(cnaqc_data_plot) = format_ID("$sID_cnaqcData")

    cnaqc_qc_plot = split_path("$cnaqc_qc_plot")
    names(cnaqc_qc_plot) = format_ID("$sID_cnaqcQC")

    table_pyclone = split_path("$pyclone_table")
    names(table_pyclone) = format_ID("$pID_table")

    pyclone_fits = split_path("$pyclone_fits")
    names(pyclone_fits) = format_ID("$pID_pycloneFit")

    pyclone_best = split_path("$pyclone_best")
    names(pyclone_best) = format_ID("$pID_pycloneBest")

    viber_pdf = split_path("$viber_pdf")
    names(viber_pdf) = format_ID("$pID_viber")

    ctree_viber_pdf = split_path("$ctree_viber_pdf")
    names(ctree_viber_pdf) = format_ID("$pID_viber")
    
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
    
    # Patient
    for (p in patient_ID){
      print(image_read_pdf(viber_pdf[[p]]) %>% image_ggplot())
      print(image_read_pdf(ctree_viber_pdf[[p]]) %>% image_ggplot())


      current_sID = sample_ID[[p]]
      # PLOT PYCLONE HARD CODED
      table =  read.table(file = table_pyclone[[p]], header = T,fill = T, sep = "\\t")
      best = read.table(pyclone_best[[p]], header = T, fill = T)
      fit = pyclone_fits[[p]]
      plot_pyclone = plot_summary_pyclone(
                      x = table,
                      y = best,
                      h5_file = fit,
                      d1 = current_sID[[1]],
                      d2 = current_sID[[2]]
                    )
      print(plot_pyclone)

      # Sample
      for (s in current_sID){
        print(image_read_pdf(cnaqc_data_plot[[s]]) %>% image_ggplot())
        print(image_read_pdf(cnaqc_qc_plot[[s]]) %>% image_ggplot())
      }    
    }
    dev.off()
    """
}
