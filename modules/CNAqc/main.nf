process CNAQC_ANALYSIS {
  publishDir params.publish_dir, mode: 'copy'

  input:
    
    tuple val(datasetID), val(patientID), val(sampleID), path(cna_RDS)
    tuple val(datasetID), val(patientID), val(sampleID), path(snv_RDS)
  
  output:

    tuple val(patientID), path("$datasetID/$patientID/$sampleID/CNAqc/*.rds"), emit: rds
    tuple val(patientID), path("$datasetID/$patientID/$sampleID/CNAqc/*.pdf"), emit: pdf
    
  script:

    def args              = task.ext.args                         ?: ''
    def matching_strategy = args!='' && args.matching_strategy    ?  "$args.matching_strategy" : "closest"

    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    library(CNAqc)
    
    res_dir = paste0("$datasetID", "/", "$patientID", "/", "$sampleID", "/CNAqc/")
    dir.create(res_dir, recursive = TRUE)

    SNV = readRDS("$snv_RDS")
    SNV = SNV[[1]]
    SNV = SNV\$mutations

    CNA = readRDS("$cna_RDS")
    x = CNAqc::init(
        mutations = SNV,
        cna = CNA\$segments,
        purity = CNA\$purity ,
        ref = "$params.assembly")
    
    # require packages not include in the singularity image
    # x = CNAqc::annotate_variants(x, drivers = CNAqc::intogen_drivers)

    x = CNAqc::analyze_peaks(x, 
      matching_strategy = "$matching_strategy")

    x = CNAqc::compute_CCF(
      x,
      muts_per_karyotype = 25,
      cutoff_QC_PASS = 0.1,
      method = "ENTROPY"
    )

    pl = ggpubr::ggarrange(
      CNAqc::plot_data_histogram(x, which = 'VAF'),
      CNAqc::plot_data_histogram(x, which = 'DP'),
      CNAqc::plot_data_histogram(x, which = 'NV'),
      CNAqc::plot_data_histogram(x, which = 'CCF'),
      ncol = 2,
      nrow = 2
    )

    pl_exp = cowplot::plot_grid(
      CNAqc::plot_gw_counts(x),
      CNAqc::plot_gw_vaf(x, N = 10000),
      CNAqc::plot_gw_depth(x, N = 10000),
      CNAqc::plot_segments(x, highlight = c("1:0", "1:1", "2:0", "2:1", '2:2')),
      pl,
      align = 'v', 
      nrow = 5,
      rel_heights = c(.15, .15, .15, .6, .5)
    ) 

    pl_qc = cowplot::plot_grid(
      CNAqc::plot_peaks_analysis(x, what = 'common', empty_plot = FALSE),
      CNAqc::plot_peaks_analysis(x, what = 'general', empty_plot = FALSE),
      CNAqc::plot_peaks_analysis(x, what = 'subclonal', empty_plot = FALSE),
      CNAqc::plot_qc(x),
      CNAqc::plot_CCF(x, assembly_plot = TRUE, empty_plot = FALSE),
      nrow = 5,
      rel_heights = c(.15, .15, .15, .3, .15)
    )

    saveRDS(object = x, file = paste0(res_dir, "qc.rds"))

    ggplot2::ggsave(plot = pl_exp, filename = paste0(res_dir, "data.pdf"), width = 12, height = 18, units = 'in', dpi = 200)
    ggplot2::ggsave(plot = pl_qc, filename = paste0(res_dir, "qc.pdf"), width = 12, height = 18, units = 'in', dpi = 200)
    """
}
