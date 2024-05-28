process VIBER {
  publishDir (
    params.publish_dir,
    mode: "copy"
  )

  input:
    
    tuple val(datasetID), val(patientID), val(sampleID), path(joint_table) // from the formatter output

  output:
    
    // tuple val(patientID), val(sampleID), path(outputPath_ctree), emit: ctree_input  // do not save or save inside mobster
    // tuple val(patientID), path("$patientID/*/mobster/*.rds")  // save also fits for each sample in mobster/sample_id/mobster.rds
    tuple val(datasetID), val(patientID), val(sampleID), path(best_fit), emit: viber_rds
    // tuple val(patientID), path("$patientID/mobster/*.pdf"), emit: mobster_pdf

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
    def n_cutoff = args!="" && args.n_cutoff ? "$args.n_cutoff" : ""
    def pi_cutoff = args!="" && args.pi_cutoff ? "$args.pi_cutoff" : ""
    def re_assign = args!="" && args.re_assign ? "$args.re_assign" : ""
    def step                    = args.step                     ?  "$args.step" : ""

      if (step == "subclonal_singlesample") {
        outDir = "$patientID/$sampleID/viber"
        outDir_ctree = "$patientID/$sampleID/ctree"
        best_fit = "$patientID/$sampleID/viber/viber_best_st_fit.rds"
      } else if (step == "subclonal_multisample"){
        sampleID = sampleID.join(' ')
        outDir = "$patientID/viber"
        outDir_ctree = "$patientID/ctree"
        best_fit = "$patientID/viber/viber_best_st_fit.rds"
        //path_ctree = "$patientID/ctree/ctree_input_pyclonevi.csv"
      }


    """
    #!/usr/bin/env Rscript

    # Sys.setenv("VROOM_CONNECTION_SIZE"=99999999)

    library(VIBER)
    library(dplyr)
    library(tidyverse)

    print("$sampleID") 
    patientID="$patientID"
    
    dir.create("$outDir", recursive = TRUE)
    samples <-strsplit(x = "$sampleID", " ")%>% unlist()
    print(samples) 
    read.csv("$joint_table", sep="\t") %>% filter(sample_id%in%samples) %>% 
      write.table(append = F,file = paste0("$outDir","/joint_table.tsv"), quote = F,sep = "\t",row.names = F)
    
    print("Subset joint done")
    ## Read input joint table
    input_tab = read.csv("$outDir/joint_table.tsv", sep="\t") %>%
      dplyr::rename(variantID = gene) %>%
      #dplyr::rename(is.driver = is_driver) %>%
      #dplyr::rename(tumour_content = purity) %>%
      dplyr::filter(patientID==patientID) %>%
      #dplyr::mutate(DP=ref_counts+alt_counts, 
      #              VAF=alt_counts/DP) %>%
      # dplyr::filter(VAF <= 1) %>% 
      dplyr::mutate(VAF=replace(VAF, VAF==0, 1e-7))
      #dplyr::rename(is_driver=is.driver) 
      #dplyr::rename(is_driver=is.driver, driver_label=variantID)

    ## Convert the input table into longer format
    reads_data = input_tab %>% 
      tidyr::pivot_wider(names_from="sample_id",
                         values_from=c("NR","NV","normal_cn",
                                       "major_cn","minor_cn","purity",
                                       "DP","VAF"), names_sep=".")

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
                                 data=reads_data,samples = 1
                                 # %>% dplyr::filter(mutation_id %in% non_tail)
                                 )

    best_fit = st_fit

    # Apply the heuristic
    st_heuristic_fit = VIBER::choose_clusters(st_fit, dimensions_cutoff = 0)
    best_fit_heuristic = st_heuristic_fit
    
    # Plot and save the fits
    if ("$step"=="subclonal_multisample"){ #mutlisample mode on
      print("multisample mode on")
      ## save fits
      saveRDS(best_fit, file = "$best_fit")
      #saveRDS(best_fit_heuristic, file = paste0("$outDir", "/viber_best_st_heuristic_fit.rds"))
      ## save plots
      plot_fit = plot(best_fit)
      plot_fit_heuristic = plot(best_fit_heuristic)
      #saveRDS(best_fit, file=paste0("$outDir", "/viber_best_st_fit.rds"))
      #saveRDS(plot_fit_heuristic, file =paste0("$outDir", "/viber_best_st_heuristic_plots.rds"))
    } else {
      print("single sample mode on")
      ## save fits
      saveRDS(best_fit, file = "$best_fit")
      #saveRDS(best_fit_heuristic, file = paste0("$outDir", "/viber_best_st_heuristic_fit.rds"))
      ## save plots
      plot_fit_mixing <- plot_mixing_proportions(best_fit)
      plot_fit_mixing_heuristic <- plot_mixing_proportions(best_fit_heuristic)
      #saveRDS(plot_fit_mixing, file = paste0("$outDir", "/viber_best_st_mixing_plots.rds"))
      #saveRDS(plot_fit_mixing_heuristic, file = paste0("$outDir", "/viber_best_st_heuristic_mixing_plots.rds"))
    }
    """
}
