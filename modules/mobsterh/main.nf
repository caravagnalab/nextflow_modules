process MOBSTERh {
  publishDir (
    params.publish_dir,
    mode: "copy"
  )

  input:
    tuple val(datasetID), val(patientID), val(sampleID), path(joint_table)

  output:
    tuple val(datasetID), val(patientID), val(sampleID), path("${outDir}/*/mobsterh_st_fit.rds"), emit: mobster_rds
    tuple val(datasetID), val(patientID), val(sampleID), path("${outDir}/*/mobsterh_st_best_fit.rds"), emit: mobster_best_rds
    tuple val(datasetID), val(patientID), val(sampleID), path("${outDir}/*/*plots.rds"), emit: mobster_plots_rds
    tuple val(datasetID), val(patientID), val(sampleID), path("${outDir}/*/REPORT_plots_mobster.rds"), emit: mobster_report_rds
    tuple val(datasetID), val(patientID), val(sampleID), path("${outDir}/*/REPORT_plots_mobster.pdf"), emit: mobster_report_pdf
    tuple val(datasetID), val(patientID), val(sampleID), path("${outDir}/*/REPORT_plots_mobster.png"), emit: mobster_report_png

  script:
    def args = task.ext.args ?: ""
    def K = args!="" && args.K ? "$args.K" : ""
    def samples = args!="" && args.samples ? "$args.samples" : ""
    def init = args!="" && args.init ? "$args.init" : ""
    def tail = args!="" && args.tail ? "$args.tail" : ""
    def epsilon = args!="" && args.epsilon ? "$args.epsilon" : ""
    def maxIter = args!="" && args.maxIter ? "$args.maxIter" : "" 
    def fit_type = args!="" && args.fit_type ? "$args.fit_type" : ""
    def seed = args!="" && args.seed ? "$args.seed" : ""
    def model_selection = args!="" && args.model_selection ? "$args.model_selection" : ""
    def trace = args!="" && args.trace ? "$args.trace" : ""
    def parallel = args!="" && args.parallel ? "$args.parallel" : ""
    def pi_cutoff = args!="" && args.pi_cutoff ? "$args.pi_cutoff" : ""
    def n_cutoff = args!="" && args.n_cutoff ? "$args.n_cutoff" : ""
    def auto_setup = args!="" && args.auto_setup ? "$args.auto_setup" : ""
    def silent = args!="" && args.silent ? "$args.silent" : ""

    outDir = "subclonal_deconvolution/mobster/$datasetID/$patientID"
    
    if (!(sampleID instanceof String)) {
      sampleID_string = sampleID.join(",")
    } else {
      sampleID_string = sampleID
    }

    """
    #!/usr/bin/env Rscript

    # Sys.setenv("VROOM_CONNECTION_SIZE"=99999999)

    library(CNAqc)
    library(mobster)
    library(dplyr)
    library(ggplot2)
    source("$moduleDir/getters.R")

    patientID = description = "$patientID"
    samples = strsplit(x = "$sampleID_string", ",") %>% unlist()  # list of samples
    dir.create("$outDir", recursive = TRUE)

    ## read mCNAqc object
    if ( grepl(".rds\$", tolower("$joint_table")) ) {
      obj = readRDS("$joint_table")
      if (class(obj) == "m_cnaqc") {
        original = obj %>% get_sample(sample=samples, which_obj="original")
        input_table = lapply(names(original), 
                             function(sample_name) 
                               CNAqc::Mutations(x=original[[sample_name]]) %>% 
                                 dplyr::mutate(sample_id=sample_name)
                             ) %>% dplyr::bind_rows()
      } else {
        cli::cli_abort("Object of class {class($joint_table)} not supported.")
      }
    } else {
      input_table = read.csv("$joint_table")
    }

    ## Function to run a single mobster fit
    run_mobster_fit = function(joint_table, descr) {
      # get input table for the single patient
      inp_tb = joint_table %>%
        dplyr::filter(VAF < 1) %>%
        dplyr::mutate(VAF=replace(VAF, VAF==0, 1e-7))

        # dplyr::rename(variantID=gene) %>%
        # dplyr::rename(is.driver=is_driver) %>%
        # dplyr::rename(tumour_content=purity) %>%
        # dplyr::rename(is_driver=is.driver) 
        # dplyr::rename(is_driver=is.driver, driver_label=variantID)

      mobster_fit(x = inp_tb,
                  K = eval(parse(text="$K")),
                  samples = as.integer("$samples"),
                  init = "$init",
                  tail = eval(parse(text="$tail")),
                  epsilon = as.numeric("$epsilon"),
                  maxIter = as.integer("$maxIter"),
                  fit.type = "$fit_type",
                  seed = as.integer("$seed"),
                  model.selection = "$model_selection",
                  trace = as.logical("$trace"),
                  parallel = as.logical("$parallel"),
                  pi_cutoff = as.numeric("$pi_cutoff"),
                  N_cutoff = as.integer("$n_cutoff"),
                  auto_setup = eval(parse(text="$auto_setup")),
                  silent = as.logical("$silent"),
                  description = descr)
    }

    lapply(samples, function(sample_name) {
      outDir_sample = paste0("$outDir/", sample_name, "/")
      dir.create(outDir_sample, recursive = TRUE)

      fit = run_mobster_fit(joint_table=input_table %>% dplyr::filter(sample_id == !!sample_name), 
                            descr=description)
      
      best_fit = fit[["best"]]
      plot_fit = plot(best_fit)

      saveRDS(object=fit, file=paste0(outDir_sample, "mobsterh_st_fit.rds"))
      saveRDS(object=best_fit, file=paste0(outDir_sample, "mobsterh_st_best_fit.rds"))
      saveRDS(object=plot_fit, file=paste0(outDir_sample, "mobsterh_st_best_fit_plots.rds"))

      # save report plots
      report_fig = mobster::plot_model_selection(fit)
      saveRDS(report_fig, file=paste0("$outDir", "REPORT_plots_mobster.rds"))
      ggplot2::ggsave(plot=report_fig, filename=paste0("$outDir", "REPORT_plots_mobster.pdf"), height=210, width=210, units="mm", dpi = 200)
      ggplot2::ggsave(plot=report_fig, filename=paste0("$outDir", "REPORT_plots_mobster.png"), height=210, width=210, units="mm", dpi = 200)
    })

    """
}
