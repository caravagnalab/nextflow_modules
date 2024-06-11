process VIBER {
  publishDir (
    params.publish_dir,
    mode: "copy"
  )

  input:
    
    tuple val(datasetID), val(patientID), val(sampleID), path(joint_table)

  output:
    
    tuple val(datasetID), val(patientID), val(sampleID), path("$outDir/viber_best_st_fit.rds"), emit: viber_rds
    tuple val(datasetID), val(patientID), val(sampleID), path("$outDir/viber_best_st_heuristic_fit.rds"), emit: viber_heuristic_rds
    tuple val(datasetID), val(patientID), val(sampleID), path("$outDir/$plot1"), emit: viber_plots_rds
    tuple val(datasetID), val(patientID), val(sampleID), path("$outDir/$plot2"), emit: viber_heuristic_plots_rds

  script:
    // viber fit params
    def args = task.ext.args ?: ""
    def K = args!="" && args.K ? "$args.K" : ""
    def alpha_0 = args!="" && args.alpha_0 ? "$args.alpha_0" : ""
    def a_0 = args!="" && args.a_0 ? "$args.a_0" : ""
    def b_0 = args!="" && args.b_0 ? "$args.b_0" : ""
    def maxIter = args!="" && args.maxIter ? "$args.maxIter" : ""
    def epsilon_conv = args!="" && args.epsilon_conv ? "$args.epsilon_conv" : ""
    def samples = args!="" && args.samples ? "$args.samples" : ""
    def q_init = args!="" && args.q_init ? "$args.q_init" : ""
    def trace = args!="" && args.trace ? "$args.trace" : ""
    // viber filter params
    def binomial_cutoff = args!="" && args.binomial_cutoff ? "$args.binomial_cutoff" : ""
    def dimensions_cutoff = args!="" && args.dimensions_cutoff ? "$args.dimensions_cutoff" : ""
    def pi_cutoff = args!="" && args.pi_cutoff ? "$args.pi_cutoff" : ""
    def re_assign = args!="" && args.re_assign ? "$args.re_assign" : ""
    def mode = args.mode ?  "$args.mode" : ""

    if (mode == "singlesample") {
      sampleID_string = sampleID
      outDir = "subclonal_deconvolution/viber/$datasetID/$patientID/$sampleID/"
      plot1 = "viber_best_st_mixing_plots.rds"
      plot2 = "viber_best_st_heuristic_mixing_plots.rds"
    } else if (mode == "multisample"){
      sampleID_string = sampleID.join(" ")
      outDir = "subclonal_deconvolution/viber/$datasetID/$patientID/"
      plot1 = "viber_best_st_fit_plots.rds"
      plot2 = "viber_best_st_heuristic_fit_plots.rds"
    }

    """
    #!/usr/bin/env Rscript

    # Sys.setenv("VROOM_CONNECTION_SIZE"=99999999)

    library(VIBER)
    library(dplyr)
    library(tidyverse)
    source("$moduleDir/getters.R")
    dir.create("$outDir", recursive = TRUE)

    patientID = "$patientID"
    samples = strsplit(x = "$sampleID_string", " ")%>% unlist()

    print("$sampleID_string")
    print("$joint_table")
    print(samples)

    if ( grepl(".rds\$", tolower("$joint_table")) ) {
      input_obj = readRDS("$joint_table")
      if (class(input_obj) == "m_cnaqc") {
        shared = input_obj %>% get_sample(sample=samples, which_obj="shared")
        joint_table = lapply(names(shared), 
                             function(sample_name) 
                               CNAqc::Mutations(x=shared[[sample_name]]) %>% 
                                 dplyr::mutate(sample_id=sample_name)
                             ) %>% dplyr::bind_rows()
      } else {
        cli::cli_alert_warning("Object of class {class(input_obj)} not supported.")
        return()
      }
    }

    print("Subset joint done")

    ## TODO : add drivers to `input_tab`

    ## Read input joint table
    input_tab = joint_table %>% 
      dplyr::mutate(VAF=replace(VAF, VAF==0, 1e-7))

    ## Convert the input table into longer format
    reads_data = input_tab %>%
      dplyr::select(chr, from, ref, alt, NV, DP, VAF, sample_id) %>% 
      tidyr::pivot_wider(names_from="sample_id",
                         values_from=c("NV","DP","VAF"), names_sep=".")

    ## Extract DP (depth)
    dp = reads_data %>%  
      # dplyr::filter(mutation_id %in% non_tail) %>% ## this step should be managed before by other module
      dplyr::select(dplyr::starts_with("DP")) %>% 
      dplyr::mutate(dplyr::across(.cols=dplyr::everything(), 
                                  .fns=function(x) replace(x, is.na(x), 0))) %>% 
      dplyr::rename_all(function(x) stringr::str_remove_all(x,"DP."))

    ## Extract NV (alt_counts)
    nv = reads_data %>% 
      # dplyr::filter(mutation_id %in% non_tail) %>% ## this step should be managed before by other module
      dplyr::select(dplyr::starts_with("NV")) %>% 
      dplyr::mutate(dplyr::across(.cols=dplyr::everything(), 
                                  .fns=function(x) replace(x, is.na(x), 0))) %>% 
      dplyr::rename_all(function(x) stringr::str_remove_all(x,"NV."))

    # Standard fit
    viber_K = eval(parse(text="$K")) 
    viber_K[which.min(viber_K)] = 2
    st_fit = VIBER::variational_fit(nv, dp, 
                                    K=viber_K, 
                                    data=reads_data,
                                    # %>% dplyr::filter(mutation_id %in% non_tail)
                                    alpha_0=as.numeric("$alpha_0"),
                                    a_0=as.integer("$a_0"),
                                    b_0=as.integer("$b_0"),
                                    max_iter=as.integer("$maxIter"),
                                    epsilon_conv=as.numeric("$epsilon_conv"),
                                    samples=as.integer("$samples"),
                                    q_init="$q_init",
                                    trace=as.logical("$trace"),
                                    description=""
                                    )

    best_fit = best_fit_heuristic = st_fit
    
    # If all clusters are removed -> keep the origianl best fit
    tryCatch(expr = {
      # Apply the heuristic
      best_fit_heuristic = VIBER::choose_clusters(st_fit, 
                                                  binomial_cutoff=as.numeric("$binomial_cutoff"),
                                                  dimensions_cutoff=as.integer("$dimensions_cutoff"),
                                                  pi_cutoff=as.numeric("$pi_cutoff"),
                                                  re_assign=as.logical("$re_assign")
                                                  )
    }, error = function(e) {
      print(e)
      best_fit_heuristic <<- st_fit
    } )

    # Save fits
    saveRDS(best_fit, file=paste0("$outDir", "viber_best_st_fit.rds"))
    saveRDS(best_fit_heuristic, file = paste0("$outDir", "viber_best_st_heuristic_fit.rds"))

    # Save plots
    if ("$mode" == "multisample") { #mutlisample mode on
      print("multisample mode on")
      plot_fit = plot(best_fit)
      plot_fit_heuristic = plot(best_fit_heuristic)
      
      saveRDS(plot_fit, file=paste0("$outDir", "viber_best_st_fit_plots.rds"))
      saveRDS(plot_fit_heuristic, file=paste0("$outDir", "viber_best_st_heuristic_fit_plots.rds"))
    } else if ("$mode" == "singlesample") {
      plot_fit_mixing = plot_mixing_proportions(best_fit)
      plot_fit_mixing_heuristic = plot_mixing_proportions(best_fit_heuristic)

      saveRDS(plot_fit_mixing, file=paste0("$outDir", "viber_best_st_mixing_plots.rds"))
      saveRDS(plot_fit_mixing_heuristic, file=paste0("$outDir", "viber_best_st_heuristic_mixing_plots.rds"))
    }
    """
}
