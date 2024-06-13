process CNAQC {
  publishDir params.publish_dir, mode: 'copy'

  input:
    
    tuple val(datasetID), val(patientID), val(sampleID), path(cna_RDS)
    tuple val(datasetID), val(patientID), val(sampleID), path(snv_RDS)
  
  output:

    tuple val(datasetID), val(patientID), val(sampleID), path("QC/CNAqc/$datasetID/$patientID/$sampleID/qc.rds"), emit: qc_rds
    tuple val(datasetID), val(patientID), val(sampleID), path("QC/CNAqc/$datasetID/$patientID/$sampleID/plot_data.rds"), path("QC/CNAqc/$datasetID/$patientID/$sampleID/plot_qc.rds"), emit: plot_rds
    tuple val(datasetID), val(patientID), val(sampleID), path("QC/CNAqc/$datasetID/$patientID/$sampleID/data.pdf"), emit: plot_pdf_data
    tuple val(datasetID), val(patientID), val(sampleID), path("QC/CNAqc/$datasetID/$patientID/$sampleID/qc.pdf"), emit: plot_pdf_qc

  script:

    def args                                = task.ext.args                                         ?: ''
    def matching_strategy                   = args!='' && args.matching_strategy                    ?  "$args.matching_strategy" : ""
    def karyotypes                          = args!='' && args.karyotypes                           ?  "$args.karyotypes" : ""
    def min_karyotype_size                  = args!='' && args.min_karyotype_size                   ?  "$args.min_karyotype_size" : ""
    def min_absolute_karyotype_mutations    = args!='' && args.min_absolute_karyotype_mutations     ?  "$args.min_absolute_karyotype_mutations" : ""
    def p_binsize_peaks                     = args!='' && args.p_binsize_peaks                      ?  "$args.p_binsize_peaks" : ""
    def matching_epsilon                    = args!='' && args.matching_epsilon                     ?  "$args.matching_epsilon" : ""
    def purity_error                        = args!='' && args.purity_error                         ?  "$args.purity_error" : ""
    def vaf_tolerance                       = args!='' && args.vaf_tolerance                        ?  "$args.vaf_tolerance" : ""
    def n_bootstrap                         = args!='' && args.n_bootstrap                          ?  "$args.n_bootstrap" : ""
    def kernel_adjust                       = args!='' && args.kernel_adjust                        ?  "$args.kernel_adjust" : ""
    def kde                                 = args!='' && args.kde                                  ?  "$args.KDE" : ""
    def starting_state_subclonal_evolution  = args!='' && args.starting_state_subclonal_evolution   ?  "$args.starting_state_subclonal_evolution" : ""
    def cluster_subclonal_CCF               = args!='' && args.cluster_subclonal_CCF                ?  "$args.cluster_subclonal_CCF" : ""

    def muts_per_karyotype                  = args!='' && args.muts_per_karyotype                   ?  "$args.muts_per_karyotype" : ""
    def cutoff_QC_PASS                      = args!='' && args.cutoff_QC_PASS                       ?  "$args.cutoff_QC_PASS" : ""
    def method                              = args!='' && args.method                               ?  "$args.method" : ""

    def plot_cn                             = args!='' && args.plot_cn                               ?  "$args.plot_cn" : ""
    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    library(CNAqc)
    
    res_dir = paste0("QC/CNAqc/", "$datasetID", "/", "$patientID", "/", "$sampleID", "/")
    dir.create(res_dir, recursive = TRUE, showWarnings = F)

    SNV = readRDS("$snv_RDS")
    SNV = SNV[["$sampleID"]]
    SNV = SNV\$mutations 
    SNV = SNV %>% dplyr::mutate(mutation_id = paste(chr,from,to,ref,alt,sep = ':'))

    CNA = readRDS("$cna_RDS")
    x = CNAqc::init(
        mutations = SNV,
        cna = CNA\$segments,
        purity = CNA\$purity ,
        ref = "$params.assembly")

    x = CNAqc::analyze_peaks(x, 
      matching_strategy = "$matching_strategy",
      karyotypes =  eval(parse(text="$karyotypes")),
      min_karyotype_size = as.numeric("$min_karyotype_size"),
      min_absolute_karyotype_mutations = as.numeric("$min_absolute_karyotype_mutations"),
      p_binsize_peaks = as.numeric("$p_binsize_peaks"),
      matching_epsilon = eval(parse(text="$matching_epsilon")),
      purity_error = as.numeric("$purity_error"),
      VAF_tolerance = as.numeric("$vaf_tolerance"),
      n_bootstrap = as.numeric("$n_bootstrap"),
      kernel_adjust = as.numeric("$kernel_adjust"),
      KDE = eval(parse(text = "$kde")),
      starting_state_subclonal_evolution = "$starting_state_subclonal_evolution",
      cluster_subclonal_CCF = as.logical("$cluster_subclonal_CCF"),
      min_VAF = 0
      )

    x = CNAqc::compute_CCF(
      x,
      karyotypes = eval(parse(text="$karyotypes")),
      muts_per_karyotype = as.numeric("$muts_per_karyotype"),
      cutoff_QC_PASS = as.numeric("$cutoff_QC_PASS"),
      method = "$method"
    )

    pl = ggpubr::ggarrange(
      CNAqc::plot_data_histogram(x, which = 'VAF', karyotypes = eval(parse(text="$karyotypes"))),
      CNAqc::plot_data_histogram(x, which = 'DP', karyotypes = eval(parse(text="$karyotypes"))),
      CNAqc::plot_data_histogram(x, which = 'NV', karyotypes = eval(parse(text="$karyotypes"))),
      CNAqc::plot_data_histogram(x, which = 'CCF', karyotypes = eval(parse(text="$karyotypes"))),
      ncol = 2,
      nrow = 2
    )

    pl_exp = cowplot::plot_grid(
      CNAqc::plot_gw_counts(x),
      CNAqc::plot_gw_vaf(x, N = 10000),
      CNAqc::plot_gw_depth(x, N = 10000),
      CNAqc::plot_segments(x, 
        highlight = eval(parse(text="$karyotypes")), 
        cn = "$plot_cn"),
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
    saveRDS(object = pl_exp, file = paste0(res_dir, "plot_data.rds"))
    saveRDS(object = pl_qc, file = paste0(res_dir, "plot_qc.rds"))

    ggplot2::ggsave(plot = pl_exp, filename = paste0(res_dir, "data.pdf"), width = 12, height = 18, units = 'in', dpi = 200)
    ggplot2::ggsave(plot = pl_qc, filename = paste0(res_dir, "qc.pdf"), width = 12, height = 18, units = 'in', dpi = 200)
    """
}
