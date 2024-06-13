process CTREE {
  publishDir params.publish_dir, mode: "copy"

  input:
    tuple val(datasetID), val(patientID), val(sampleID), path(ctree_input)

  output:
    tuple val(datasetID), val(patientID), val(sampleID), path("${outDir}/ctree_{VIBER,MOBSTERh,pyclonevi}.rds"), emit: ctree_rds, optional: true 
    tuple val(datasetID), val(patientID), val(sampleID), path("${outDir}/ctree_{VIBER,MOBSTERh,pyclonevi}_plots.rds"), emit: ctree_plots_rds, optional: true
    tuple val(datasetID), val(patientID), val(sampleID), path("${outDir}/REPORT_plots_ctree_{VIBER,MOBSTERh,pyclonevi}.rds"), emit: ctree_report_rds, optional: true
    tuple val(datasetID), val(patientID), val(sampleID), path("${outDir}/REPORT_plots_ctree_{VIBER,MOBSTERh,pyclonevi}.pdf"), emit: ctree_report_pdf, optional: true
    tuple val(datasetID), val(patientID), val(sampleID), path("${outDir}/REPORT_plots_ctree_{VIBER,MOBSTERh,pyclonevi}.png"), emit: ctree_report_png, optional: true

  script:

    def args = task.ext.args ?: ""
    def sspace_cutoff = args!="" && args.sspace_cutoff ? "$args.sspace_cutoff" : ""
    def n_sampling = args!="" && args.n_sampling ? "$args.n_sampling" : ""
    def store_max = args!="" && args.store_max ? "$args.store_max" : ""
    def mode = args!="" && args.mode ?  "$args.mode" : ""

    if (mode == "singlesample") {
      outDir = "subclonal_deconvolution/ctree/$datasetID/$patientID/$sampleID"
    } else if (mode == "multisample") {
      outDir = "subclonal_deconvolution/ctree/$datasetID/$patientID"
    }

    """
    #!/usr/bin/env Rscript

    library(ctree)
    library(dplyr)
    library(VIBER)
    library(mobster)
    
    initialize_ctree_obj = function(ctree_input) {
      if (!"variantID" %in% colnames(ctree_input) | !"is.driver" %in% colnames(ctree_input)) {
        ctree_input = ctree_input %>% dplyr::mutate(is.driver=FALSE, variantID=NA)
        ctree_input[1, "is.driver"] = TRUE; ctree_input[1, "variantID"] = ""
      }

      driver_cluster = unique(ctree_input[which(ctree_input["is.driver"]==TRUE),c("cluster")])
      # the CCF table must report CCF values for each cluster and sample
      # cluster | nMuts | is.driver | is.clonal | sample1 | sample2 | ...
      CCF_table = ctree_input %>% 
        dplyr::select(sample_id, cluster, nMuts, is.driver, is.clonal, CCF) %>%
        dplyr::filter(is.driver  != "") %>%
        dplyr::filter(cluster != "Tail") %>%
        mutate(cluster = as.character(cluster)) %>%
        group_by(cluster) %>%
        dplyr::filter(!(cluster %in% driver_cluster) | is.driver  == "TRUE") %>%
        mutate(is.driver = as.logical(is.driver)) %>%
        dplyr::ungroup() %>%
        unique() %>%
        tidyr::pivot_wider(names_from = "sample_id", values_from = "CCF", values_fill = 0)
      
      # the driver table must contain patient and variant IDs and report clonality and driver status
      # patientID | variantID | is.driver | is.clonal | cluster | sample1 | sample2 | ...
      drivers_table = ctree_input %>% 
         dplyr::select(patientID, sample_id, variantID, cluster, is.driver, is.clonal, CCF) %>%
         dplyr::filter(is.driver==TRUE) %>%
         mutate(is.driver = as.logical(is.driver)) %>%
         dplyr::mutate(cluster = as.character(cluster)) %>%
         tidyr::pivot_wider(names_from="sample_id", values_from="CCF",values_fill = 0)

      samples = unique(ctree_input[["sample_id"]])  # if multisample, this is a list
      patient = unique(ctree_input[["patientID"]])
      ctree_init = list("CCF_table"=CCF_table,
                        "drivers_table"=drivers_table,
                        "samples"=samples,
                        "patient"=patient)
      return(ctree_init)
    }

    if ( grepl(".rds\$", tolower("$ctree_input")) ) {
      best_fit = readRDS("$ctree_input")
      do_fit = TRUE

      if (!"data" %in% names(best_fit)) { best_fit[["data"]] =  data.frame(row.names=1:best_fit[["N"]])}

      ## viber
      if (class(best_fit) == "vb_bmm") {
        fn_name = VIBER::get_clone_trees
        subclonal_tool = "VIBER"
        
        # If only one sample inference ctree is not working with VIBER
        if (ncol(best_fit[["x"]]) <= 2) {
          cli::cli_alert_warning("Object of class {class(best_fit)} with only one sample is not supported by the function `get_clone_trees`!")
          do_fit = FALSE
        }

        if (!"gene" %in% colnames(best_fit[["data"]]) | !"driver" %in% colnames(best_fit[["data"]])) {
          best_fit[["data"]] = best_fit[["data"]] %>% dplyr::mutate(driver=FALSE, gene=NA)
          best_fit[["data"]][1, "driver"] = TRUE; best_fit[["data"]][1, "gene"] = ""
        }
      }

      ## mobster
      if (class(best_fit) == "dbpmm") {
        fn_name = mobster::get_clone_trees
        subclonal_tool = "MOBSTERh"
        if (!"driver_label" %in% colnames(best_fit[["data"]]) | !"is_driver" %in% colnames(best_fit[["data"]])) {
          best_fit[["data"]] = best_fit[["data"]] %>% dplyr::mutate(is_driver=FALSE, driver_label=NA)
          best_fit[["data"]][1, "is_driver"] = TRUE; best_fit[["data"]][1, "driver_label"] = ""
        }
      }

      if (class(best_fit) %in% c("vb_bmm", "dbpmm") & do_fit) {
        # VIBER or MOBSTER object
        trees = fn_name(x = best_fit,
                        sspace.cutoff = as.integer("$sspace_cutoff"),
                        n.sampling = as.integer("$n_sampling"),
                        store.max = as.integer("$store_max"))
      } else {
        cli::cli_alert_warning("Object of class {class(best_fit)} not supported.")
        do_fit = FALSE
      }

    } else {
      do_fit = TRUE
      subclonal_tool = "pyclonevi"
      input_table = read.csv("$ctree_input", sep="\t")
      data_ctree = initialize_ctree_obj(input_table)
      trees = ctrees(CCF_clusters = data_ctree[["CCF_table"]],
                     drivers = data_ctree[["drivers_table"]],
                     samples = data_ctree[["samples"]],
                     patient = data_ctree[["patient"]],
                     sspace.cutoff = as.integer("$sspace_cutoff"),
                     n.sampling = as.integer("$n_sampling"),
                     store.max = as.integer("$store_max"))
    }

    if (do_fit) {
      ctree_output = paste0("ctree_", subclonal_tool)

      # plot the best tree
      plot_tree = plot(trees[[1]])

      # save rds and plots
      dir.create("$outDir/", recursive = TRUE)
      saveRDS(object=trees, file=paste0("$outDir/", ctree_output, ".rds"))
      saveRDS(object=plot_tree, file=paste0("$outDir/", ctree_output, "_plots.rds"))

      # Save report plot
      #top_phylo = plot(trees[[1]])
      #phylos = ggplot()
      #if (length(trees) > 1) phylos = lapply(trees[2:min(length(trees), 3)], plot) %>% patchwork::wrap_plots(nrow=1)
      #ccf = plot_CCF_clusters(trees[[1]])
      #info_transfer = plot_information_transfer(trees[[1]])
      #clone_size = plot_clone_size(trees[[1]])

      #report_fig = patchwork::wrap_plots(ccf, info_transfer, top_phylo, clone_size, phylos, design="A#BB\nD#CC\nEEEE")
      #saveRDS(object=report_fig, file=paste0("$outDir/REPORT_plots_", ctree_output, ".rds"))
      #ggsave(plot=report_fig, filename=paste0("$outDir/REPORT_plots_", ctree_output, ".pdf"), height=297, width=210, units="mm")
      #ggsave(plot=report_fig, filename=paste0("$outDir/REPORT_plots_", ctree_output, ".png"), height=297, width=210, units="mm")
    }

    """
}
