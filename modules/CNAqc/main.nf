process CNAQC {
  publishDir params.publish_dir, mode: 'copy'

  input:
    
    tuple val(datasetID), val(patientID), val(sampleID), path(cna_RDS)
    tuple val(datasetID), val(patientID), val(sampleID), path(snv_RDS)
  
  output:

    tuple val(datasetID), val(patientID), val(sampleID), path("QC/CNAqc/$datasetID/$patientID/$sampleID/*.rds"), emit: rds
    tuple val(datasetID), val(patientID), val(sampleID), path("QC/CNAqc/$datasetID/$patientID/$sampleID/*.pdf"), emit: pdf
    
  script:

    def args                                = task.ext.args                                         ?: ''
    def matching_strategy                   = args!='' && args.matching_strategy                    ?  "$args.matching_strategy" : "closest"
    def karyotypes                          = args!='' && args.karyotypes                           ?  "$args.karyotypes" : "c("1:0", "1:1", "2:0", "2:1", "2:2")"
    def min_karyotype_size                  = args!='' && args.min_karyotype_size                   ?  "$args.min_karyotype_size" : "0"
    def min_absolute_karyotype_mutations    = args!='' && args.min_absolute_karyotype_mutations     ?  "$args.min_absolute_karyotype_mutations" : "100"
    def p_binsize_peaks                     = args!='' && args.p_binsize_peaks                      ?  "$args.p_binsize_peaks" : "0.005"
    def matching_epsilon                    = args!='' && args.matching_epsilon                     ?  "$args.matching_epsilon" : "NULL"
    def purity_error                        = args!='' && args.purity_error                         ?  "$args.purity_error" : "0.05"
    def VAF_tolerance                       = args!='' && args.VAF_tolerance                        ?  "$args.VAF_tolerance" : "0.015"
    def n_bootstrap                         = args!='' && args.n_bootstrap                          ?  "$args.n_bootstrap" : "1"
    def kernel_adjust                       = args!='' && args.kernel_adjust                        ?  "$args.kernel_adjust" : "1"
    def kde                                 = args!='' && args.kde                                  ?  "$args.KDE" : "TRUE"
    def starting_state_subclonal_evolution  = args!='' && args.starting_state_subclonal_evolution   ?  "$args.starting_state_subclonal_evolution" : "1:1"
    def cluster_subclonal_CCF               = args!='' && args.cluster_subclonal_CCF                ?  "$args.cluster_subclonal_CCF" : "FALSE"

    def muts_per_karyotype                  = args!='' && args.muts_per_karyotype                   ?  "$args.muts_per_karyotype" : "25"
    def cutoff_QC_PASS                      = args!='' && args.cutoff_QC_PASS                       ?  "$args.cutoff_QC_PASS" : "0.1"
    def method                              = args!='' && args.method                               ?  "$args.method" : "ENTROPY"

    def plot_cn                             = args!='' && args.plot_cn                               ?  "$args.plot_cn" : "absolute"
    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    library(CNAqc)
    
    res_dir = paste0("QC/CNAqc/", "$datasetID", "/", "$patientID", "/", "$sampleID", "/")
    dir.create(res_dir, recursive = TRUE, showWarnings = F)

    SNV = readRDS("$snv_RDS")
    SNV = SNV[["$sampleID"]]
    SNV = SNV\$mutations

    CNA = readRDS("$cna_RDS")
    x = CNAqc::init(
        mutations = SNV,
        cna = CNA\$segments,
        purity = CNA\$purity ,
        ref = "$params.assembly")
    
    x\$mutations = x\$mutations %>% filter(VAF > 0)

    x = CNAqc::analyze_peaks(x, 
      matching_strategy = "$matching_strategy",
      karyotypes =  eval(parse(text="$karyotypes")),
      min_karyotype_size = as.integer("$min_karyotype_size"),
      min_absolute_karyotype_mutations = as.integer("$min_absolute_karyotype_mutations"),
      p_binsize_peaks = as.integer("$p_binsize_peaks"),
      matching_epsilon = eval(parse(text="$matching_epsilon")),
      purity_error = as.integer("$purity_error"),
      VAF_tolerance = as.integer("$VAF_tolerance"),
      n_bootstrap = as.integer("$n_bootstrap"),
      kernel_adjust = as.integer("$kernel_adjust"),
      KDE = as.logical("$kde"),
      starting_state_subclonal_evolution = "$starting_state_subclonal_evolution",
      cluster_subclonal_CCF = as.logical( "$cluster_subclonal_CCF"),
      #min_VAF = 0,
      )

    x = CNAqc::compute_CCF(
      x,
      karyotypes = eval(parse(text="$karyotypes")),
      muts_per_karyotype = as.integer("$muts_per_karyotype"),
      cutoff_QC_PASS = as.integer("$cutoff_QC_PASS"),
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

    ggplot2::ggsave(plot = pl_exp, filename = paste0(res_dir, "data.pdf"), width = 12, height = 18, units = 'in', dpi = 200)
    ggplot2::ggsave(plot = pl_qc, filename = paste0(res_dir, "qc.pdf"), width = 12, height = 18, units = 'in', dpi = 200)
    """
}
