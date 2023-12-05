process VIBER {

  publishDir params.publish_dir, mode: 'copy'
  
  input:
    
    tuple val(patientID), path(joint_table)
  
  output:
  
    tuple path("$patientID/viber/*pdf"), path("$patientID/viber/best_fit.rds")
  
  script:
    def args = task.ext.args ?: ''
    def K = args!='' && args.K ? "$args.K" : ""
    def samples = args!='' && args.samples ? "$args.samples" : ""
    def alpha_0 = args!='' && args.alpha_0 ? "$args.alpha_0" : ""
    def a_0 = args!='' && args.a_0 ? "$args.a_0" : ""
    def b_0 = args!='' && args.b_0 ? "$args.b_0" : ""
    def maxIter = args!='' && args.maxIter ? "$args.maxIter" : ""
    def epsilon_conv = args!='' && args.epsilon_conv ? "$args.epsilon_conv" : ""
    def q_init = args!='' && args.q_init ? "$args.q_init" : ""
    def trace = args!='' && args.trace ? "$args.trace" : ""
    def binomial_cutoff = args!='' && args.binomial_cutoff ? "$args.binomial_cutoff" : ""
    def dimensions_cutoff = args!='' && args.dimensions_cutoff ? "$args.dimensions_cutoff" : ""
    def pi_cutoff = args!='' && args.pi_cutoff ? "$args.pi_cutoff" : ""
    def re_assign = args!='' && args.re_assign ? "$args.re_assign" : ""

    """
    #!/usr/bin/env Rscript

    Sys.setenv("VROOM_CONNECTION_SIZE"=99999999)

    library(tidyverse)
    library(VIBER)
    library(dplyr)
    dir.create(paste0("$patientID","/","viber"), recursive = TRUE)

    pyclone_tsv = read.table(file = "$joint_table", header = T)
    pyclone_tsv["tot_counts"] <- pyclone_tsv["ref_counts"]+pyclone_tsv["alt_counts"]
    mobster_input = list(
     "successes"=pyclone_tsv %>% dplyr::select(mutation_id, sample_id, alt_counts) %>%
     tidyr::pivot_wider(id_cols="mutation_id", names_from="sample_id", values_from="alt_counts", values_fill=0) %>%
     tibble::column_to_rownames(var="mutation_id"),
     "trials"=pyclone_tsv %>% dplyr::select(mutation_id, sample_id, tot_counts) %>%
     tidyr::pivot_wider(id_cols="mutation_id", names_from="sample_id", values_from="tot_counts", values_fill=0) %>%
     tibble::column_to_rownames(var="mutation_id")
    )
    fit = variational_fit(x = mobster_input[["successes"]], 
                          y = mobster_input[["trials"]],
                          data = NULL,
                          K = as.integer("$K"),
                          alpha_0 = as.numeric("$alpha_0"),
                          a_0 = as.numeric("$a_0"),
                          b_0 = as.numeric("$b_0"),
                          max_iter = as.integer("$maxIter"),
                          epsilon_conv = as.numeric("$epsilon_conv"),
                          samples = as.integer("$samples"),
                          q_init = "$q_init",
                          trace = as.logical("$trace"),
                          description = "$patientID")
                          
    fit = choose_clusters(fit, 
                          binomial_cutoff = as.numeric("$binomial_cutoff"), 
                          dimensions_cutoff = as.integer("$dimensions_cutoff"),
                          pi_cutoff = as.numeric("$pi_cutoff"),
                          re_assign = as.logical("$re_assign"))

    saveRDS(object = fit, file = paste0("$patientID","/","viber","/","best_fit.rds"))
    plot_fit <- plot(fit)
    plot_mixing <- plot_mixing_proportions(fit)
    #plot_LV < plot_latent_variables(fit)
    #plot_elbo <-  plot_ELBO(fit)
    #plot_peaks <- plot_peaks(fit)
    
    pdf(paste0("$patientID","/","viber","/","plot_fit.pdf"), height=6, width=6)
    print(plot_fit)
    print(plot_mixing)
    #print(plot_LV)
    #print(plot_elbo)
    #print(plot_peaks)
    dev.off()

    """
}
