process CTREE {
  publishDir params.publish_dir, mode: "copy"

  input:
    tuple val(datasetID), val(patientID), val(sampleID), path(ctree_input)
    // val(subclonal_tool) // tool used for subclonal deconvolution

  output:
    tuple val(datasetID), val(patientID), val(sampleID), path("$patientID/ctree/*.rds"), emit: ctree_rds 
    tuple val(datasetID), val(patientID), val(sampleID), path("$patientID/ctree/*.pdf"), emit: ctree_plot 

  script:
    def args = task.ext.args ?: ""
    def sspace_cutoff = args!="" && args.sspace_cutoff ? "$args.sspace_cutoff" : ""
    def n_sampling = args!="" && args.n_sampling ? "$args.n_sampling" : ""
    def store_max = args!="" && args.store_max ? "$args.store_max" : ""
    
    if (mode == "singlesample") {
      outDir = "subclonal_deconvolution/ctree/$datasetID/$patientID/$sampleID/"
    } else if (mode == "multisample"){
      sampleID = sampleID.join(" ")
      outDir = "subclonal_deconvolution/ctree/$datasetID/$patientID/"
    }

    """
    #!/usr/bin/env Rscript

    library(ctree)
    library(dplyr)
    library(viber)
    library(mobster)
    
    initialize_ctree_obj = function(ctree_input){
      # the CCF table must report CCF values for each cluster and sample
      # cluster | nMuts | is.driver | is.clonal | sample1 | sample2 | ...
      CCF_table = ctree_input %>% 
        dplyr::select(sampleID, cluster, nMuts, is.driver, is.clonal, CCF) %>% 
        dplyr::filter(cluster != "Tail") %>% 
        mutate(cluster = as.character(cluster)) %>%
        unique() %>% 
        tidyr::pivot_wider(names_from="sampleID", values_from="CCF",values_fill = 0)
      
      # the driver table must contain patient and variant IDs and report clonality and driver status
      # patientID | variantID | is.driver | is.clonal | cluster | sample1 | sample2 | ...
      drivers_table = ctree_input %>% 
        dplyr::select(patientID, sampleID, variantID, cluster, is.driver, is.clonal, CCF) %>% 
        dplyr::filter(is.driver) %>% 
        mutate(cluster = as.character(cluster)) %>%
        tidyr::pivot_wider(names_from="sampleID", values_from="CCF",values_fill = 0)
      
      samples = unique(ctree_input[["sampleID"]])  # if multisample, this is a list
      patient = unique(ctree_input[["patientID"]])
      ctree_init = list("CCF_table"=CCF_table,
                        "drivers_table"=drivers_table,
                        "samples"=samples,
                        "patient"=patient)
      return(ctree_init)
    }

    if (tolower(file_ext($ctree_input)) == "rds") {
      best_fit = readRDS($ctree_input)
      if (class(best_fit) %in% c("vb_bmm", "dbpmm")) {
        # VIBER or MOBSTER object
        trees = get_clone_trees(x = best_fit,
                                sspace.cutoff = as.integer("$sspace_cutoff"),
                                n.sampling = as.integer("$n_sampling"),
                                store.max = as.integer("$store_max"))
        ctree_output = dplyr::case_when(
          class(best_fit) == "vb_bmm" ~ subclonal_tool = "VIBER",
          class(best_fit) == "dbpmm" ~ subclonal_tool = "MOBSTERh"
        )
      } else {
        cli::cli_alert_warning("Object of class {class($ctree_input)} not supported.")
        return()
      }
    } else {
      subclonal_tool = "pyclonevi"
      input_table = read.csv($ctree_input, sep="\t")
      data_ctree = initialize_ctree_obj(input_table)
      trees = ctrees(CCF_clusters = data_ctree[["CCF_table"]],
                     drivers = data_ctree[["drivers_table"]],
                     samples = data_ctree[["samples"]],
                     patient = data_ctree[["patient"]],
                     sspace.cutoff = as.integer("$sspace_cutoff"),
                     n.sampling = as.integer("$n_sampling"),
                     store.max = as.integer("$store_max"))
    }

    ctree_output = paste0("ctree_", subclonal_tool, ".rds")

    # plot the best tree
    plot_tree = plot(trees[[1]]) 

    # save rds and plots
    dir.create($outDir, recursive = TRUE)
    saveRDS(object=trees, file=paste0($outDir, ctree_output))
    """
}
