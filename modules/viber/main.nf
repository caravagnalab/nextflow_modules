process VIBER {

  publishDir params.publish_dir, mode: 'copy'
  
  input:
    
    tuple val(patientID), path(joint_table)
  
  output:
  
    tuple path("$patientID/viber/*pdf"), path("$patientID/viber/best_fit.rds")
  
  script:
    
    """
    #!/usr/bin/env Rscript

    Sys.setenv("VROOM_CONNECTION_SIZE"=99999999)
    
    library(tidyverse)
    library(VIBER)
    dir.create(paste0("$patientID","/","viber"), recursive = TRUE)

    pyclone_tsv = read.table(file = "$joint_table", header = T)
    pyclone_tsv["tot_counts"] <- pyclone_tsv["ref_counts"]+pyclone_tsv["alt_counts"]
    mobster_input <- list()
    n_samples = length(unique(pyclone_tsv["sample_id"]))
    mobster_input_success<- list()
    mobster_input_trials<- list()

    for (s in seq(n_samples)){
        sampleID = unique(pyclone_tsv["sample_id"])[s]
        n_success_0 = pyclone_tsv%>%filter(sample_id==sampleID)%>%select(mutation_id,alt_counts)
        n_success = n_success_0["alt_counts"]
        names(n_success) <- n_success_0["mutation_id"]
        n_trials_0 = pyclone_tsv%>%filter(sample_id==sampleID)%>%select(mutation_id,tot_counts)
        n_trials = n_trials_0["tot_counts"]
        names(n_trials ) <- n_trials_0["mutation_id"]

        mobster_input_success[[s]]<-n_success
        mobster_input_trials[[s]]<-n_trials
    }
    mobster_succes <-as.data.frame(do.call(cbind, mobster_input_success))
    mobster_trials <-as.data.frame(do.call(cbind, mobster_input_trials))
    mobster_input <- list(mobster_succes, mobster_trials)
    names(mobster_input)<- c("successes","trials")

    fit = variational_fit( mobster_input[["successes"]]%>%select(V1,V2),
        mobster_input[["trials"]]%>%select(V1,V2)
        )
    fit = choose_clusters(fit, 
                      binomial_cutoff = 0, 
                      dimensions_cutoff = 0,
                      pi_cutoff = 0.02)
    saveRDS(object = fit, file = paste0("$patientID","/","viber","/","best_fit.rds"))
    plot_fit <- plot(fit)
    plot_mixing <- plot_mixing_proportions(fit)
    plot_LV < plot_latent_variables(fit)
    plot_elbo <-  plot_ELBO(fit)
    plot_peaks <- plot_peaks(fit)
    ggsave(filename = paste0("$patientID","/","viber","/","plot_fit.pdf"), plot = plot_fit)
    """
}
