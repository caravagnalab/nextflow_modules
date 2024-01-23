process SEQUENZA_CNAqc {

  publishDir params.publish_dir, mode: 'copy'
  
  input:
    
    tuple val(datasetID), val(patientID), val(sampleID), val(sex), path(seqzFile), path(snv_vcfFile)
  
  output:
  
    tuple(path("$patientID/$sampleID/SEQUENZA_CNAqc/final"), path("$patientID/$sampleID/SEQUENZA_CNAqc/*.rds"), path("$patientID/$sampleID/SEQUENZA_CNAqc/*.pdf"))
  
  script:
  
    def args              = task.ext.args                         ?: ''
    def norm_method       = args!='' && args.norm_method          ? "$args.norm_method" : "median"
    def window            = args!='' && args.window               ? "$args.window" : "1e5" 
    def gamma             = args!='' && args.gamma                ? "$args.gamma" : "280" 
    def kmin              = args!='' && args.kmin                 ? "$args.kmin" : "300"  
    def min_reads_baf     = args!='' && args.min_reads_baf        ? "$args.min_reads_baf" : "50"
    def min_reads         = args!='' && args.min_reads            ? "$args.min_reads" : "50"
    def min_reads_normal  = args!='' && args.min_reads_normal     ? "$args.min_reads_normal": "15"
    def max_mut_types     = args!='' && args.max_mut_types        ? "$args.max_mut_types" : "1"
    
    def low_cell          = args!='' && args.low_cell             ? "$args.low_cell" : "0.9"
    def up_cell           = args!='' && args.up_cell              ? "$args.up_cell" : "1.0"
    def low_ploidy        = args!='' && args.low_ploidy           ? "$args.low_ploidy" : "1.8"
    def up_ploidy         = args!='' && args.up_ploidy            ? "$args.up_ploidy" : "5.4"
    def delta_cellularity = args!='' && args.delta_cellularity    ? "$args.delta_cellularity" : "0.05"
    def delta_ploidy      = args!='' && args.delta_ploidy         ? "$args.delta_ploidy" : "0.25"

    def matching_strategy = args!='' && args.matching_strategy    ? "$args.matching_strategy" : "closest"

    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    library(sequenza)
    library(CNAqc)
    library(evoverse)
    
    res_dir = paste0("$datasetID", "/", "$patientID", "/", "$sampleID", "/SEQUENZA_CNAqc/")
    dir.create(res_dir, recursive = TRUE)

    SNV = evoverse::evoparse_mutect_mutations("$snv_vcfFile")
    SNV = SNV[[1]] 
    SNV\$mutations = SNV\$mutations %>% dplyr::select(chr, from, to, ref, alt, NV, DP, VAF)

    file = R.utils::gunzip("$seqzFile", remove=FALSE)
    out = CNAqc::Sequenza_CNAqc(
      sample_id = paste0("$patientID"),
      seqz_file = file,
      mutations = SNV\$mutations,
      sex = "$sex",
      cellularity = c(as.integer("$low_cell"), as.integer("$up_cell")),
      ploidy = c(as.integer("$low_ploidy"), as.integer("$up_ploidy")),
      reference = "$params.assembly",
      normalization.method = "$norm_method",
      window = as.numeric("$window"),
      gamma = as.integer("$gamma"),
      kmin = as.integer("$kmin"),
      min.reads.baf = as.integer("$min_reads_baf"),
      min.reads = as.integer("$min_reads"),
      min.reads.normal = as.integer("$min_reads_normal"),
      max.mut.types = as.integer("$max_mut_types"),
      delta_cellularity = as.numeric("$delta_cellularity"),
      delta_ploidy = as.numeric("$delta_ploidy"),
      matching_strategy = "$matching_strategy"
    )
    
    best_solution = out %>% dplyr::filter(QC == 'PASS') %>% dplyr::filter(score == min(score)) 
    x = best_solution\$cnaqc[[1]] 
    x = CNAqc::compute_CCF(
      x,
      muts_per_karyotype = 25,
      cutoff_QC_PASS = 0.1,
      method = "ENTROPY"
    )

    best_solution\$cnaqc[[1]] = x 

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

    file.copy(paste0("./final"), res_dir, recursive=TRUE)
    saveRDS(object = best_solution, file = paste0(res_dir, "best_solution.rds"))

    ggplot2::ggsave(plot = pl_exp, filename = paste0(res_dir, "data.pdf"), width = 12, height = 18, units = 'in', dpi = 200)
    ggplot2::ggsave(plot = pl_qc, filename = paste0(res_dir, "qc.pdf"), width = 12, height = 18, units = 'in', dpi = 200)

    """
}
