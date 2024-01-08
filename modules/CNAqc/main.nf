process CNAqc {
  publishDir params.publish_dir, mode: 'copy'

  input:
    
    tuple val(datasetID), val(patientID), val(sampleID), path(cnaDir), path(snv_vcfFile)
  
  output:
  
    tuple(path("$datasetID/$patientID/$sampleID/CNAqc/*.rds"), path("$datasetID/$patientID/$sampleID/CNAqc/*.pdf"))
  
  script:

    def args              = task.ext.args                         ?: ''
    def cna_caller        = args!='' && args.cna_caller           ? "$args.cna_caller" : "sequenza"
    def variant_caller    = args!='' && args.variant_caller       ? "$args.variant_caller" : "mutect"
    def matching_strategy = args!='' && args.matching_strategy    ? "$args.matching_strategy" : "rightmost"

    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    library(sequenza)
    library(CNAqc)
    library(evoverse)
    
    res_dir = paste0("$datasetID", "/", "$patientID", "/", "$sampleID", "/CNAqc/")
    dir.create(res_dir, recursive = TRUE)

    if ("$variant_caller" == "mutect"){
        SNV = evoverse::evoparse_mutect_mutations("$snv_vcfFile")
        SNV = SNV[[1]]
        SNV = SNV\$mutations

    } else if ("$variant_caller" == "platypus"){
        SNV = evoverse::evoparse_platypus_mutations("$snv_vcfFile")
        SNV = SNV[[1]]
        SNV = SNV\$mutations
    }

    SNV = SNV %>% 
        dplyr::filter(ref %in% c('A', 'C', 'T', 'G'), alt %in% c('A', 'C', 'T', 'G')) %>%
        dplyr::filter(VAF > 0.05) %>% 
        dplyr::select(chr, from, to, ref, alt, NV, DP, VAF)

    if ("$cna_caller" == "sequenza"){
        CNA = evoverse::evoparse_Sequenza_CNAs(paste("$cnaDir", "$patientID", sep = '/'))
    }

    x = CNAqc::init(
        mutations = SNV,
        cna = CNA\$segments,
        purity = CNA\$purity ,
        ref = "$params.assembly")

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
