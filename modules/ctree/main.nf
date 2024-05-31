process CTREE {
  publishDir params.publish_dir, mode: 'copy'

  input:
    tuple val(datasetID), val(patientID), val(sampleID), path(ctree_input)
    val(subclonal_tool) // tool used for subclonal deconvolution

  output:
    tuple val(datasetID), val(patientID), val(sampleID), path("$patientID/ctree/*.rds"), emit: ctree_rds 
    tuple val(datasetID), val(patientID), val(sampleID), path("$patientID/ctree/*.pdf"), emit: ctree_plot 

  script:
    def args = task.ext.args ?: ''
    def sspace_cutoff = args!='' && args.sspace_cutoff ? "$args.sspace_cutoff" : ""
    def n_sampling = args!='' && args.n_sampling ? "$args.n_sampling" : ""
    def store_max = args!='' && args.store_max ? "$args.store_max" : ""
    
    if (mode == "singlesample") {
        outDir = "subclonal_deconvolution/ctree/$datasetID/$patientID/$sampleID/"
        // output_ctree = "subclonal_deconvolution/ctree/$datasetID/$patientID/$sampleID/ctree.rds"
      } else if (mode == "multisample"){
        sampleID = sampleID.join(' ')
        outDir = "subclonal_deconvolution/ctree/$datasetID/$patientID/"
        // output_ctree = "subclonal_deconvolution/ctree/$datasetID/$patientID/ctree.rds"
      }
    """
    #!/usr/bin/env Rscript

    out_dirname = $outDir

    library(ctree)
    library(dplyr)
    library(viber)
    library(mobster)
    
    initialize_ctree_obj <- function(ctree_input){
  
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
      ctree_init <- list("CCF_table"=CCF_table,
                        "drivers_table"=drivers_table,
                        "samples"=samples,
                        "patient"=patient)
      return(ctree_init)
    }

    if ($subclonal_tool == "pyclonevi"){
      ctree_output = paste0("ctree_",$subclonal_tool,".rds")
      data_ctree <- initialize_ctree_obj($ctree_input)
      trees = ctrees(CCF_clusters = data_ctree[["CCF_table"]],
                    drivers = data_ctree[["drivers_table"]],
                    samples = data_ctree[["samples"]],
                    patient = data_ctree[["patient"]],
                    sspace.cutoff = as.integer("$sspace_cutoff"),
                    n.sampling = as.integer("$n_sampling"),
                    store.max = as.integer("$store_max"))
    } else {
      ctree_output = paste0("ctree_",$subclonal_tool,".rds")
      best_fit = readRDS($ctree_input)
      trees = get_clone_trees(x = best_fit,
                              sspace.cutoff = as.integer("$sspace_cutoff"),
                              n.sampling = as.integer("$n_sampling"),
                              store.max = as.integer("$store_max"))
    }
    # plot the best tree
    plot_tree = plot(trees[[1]]) 

    # save rds and plots
    dir.create(out_dirname, recursive = TRUE)
    
    saveRDS(object=trees, file=paste0(out_dirname, ctree_output))






    #################################################### OLD SCRIPT #################################################
    ctree_input = read.csv("$ctree_input")
    idd = ""
    if ("tool" %in% colnames(ctree_input)) idd = paste0("_", unique(ctree_input[["tool"]]))

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

    trees = ctrees(CCF_clusters = CCF_table,
                  drivers = drivers_table,
                  samples = samples,
                  patient = patient,
                  sspace.cutoff = as.integer("$sspace_cutoff"),
                  n.sampling = as.integer("$n_sampling"),
                  store.max = as.integer("$store_max"))

    # plot the best tree
    plot_tree = plot(trees[[1]]) 

    # save rds and plots
    dir.create(out_dirname, recursive = TRUE)
    
    saveRDS(object=trees, file=paste0(out_dirname, "ctree_obj", idd, ".rds"))

    #pdf(paste0(out_dirname, "ctree_plot", idd, ".pdf"))
    #print(plot_tree)
    #dev.off()
    """
}
