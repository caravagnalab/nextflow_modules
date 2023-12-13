process CTREE {
  publishDir params.publish_dir, mode: 'copy'

  input:
    tuple val(patientID), val(timepointID), val(sampleID), path(ctree_input)

  output:
    tuple path("$patientID/$timepointID/$sampleID/*.rds"), path("$patientID/$timepointID/$sampleID/*.pdf")

  script:
    def args = task.ext.args ?: ''
    def sspace_cutoff = args!='' && args.sspace_cutoff ? "$args.sspace_cutoff" : ""
    def n_sampling = args!='' && args.n_sampling ? "$args.n_sampling" : ""
    def store_max = args!='' && args.store_max ? "$args.store_max" : ""

    """
    #!/usr/bin/env Rscript

    out_dirname = paste0("$patientID","/","$timepointID","/","$sampleID", "/")

    library(ctree)
    library(dplyr)

    ctree_input = read.csv("$ctree_input")

    # the CCF table must report CCF values for each cluster and sample
    # cluster | nMuts | is.driver | is.clonal | sample1 | sample2 | ...
    CCF_table = ctree_input %>% 
      dplyr::filter(cluster != "Tail") %>% 
      dplyr::select(-variantID) %>% unique() %>% 
      tidyr::pivot_wider(names_from="sampleID", values_from="CCF")

    # the driver table must contain patient and variant IDs and report clonality and driver status
    # patientID | variantID | is.driver | is.clonal | cluster | sample1 | sample2 | ...
    drivers_table = ctree_input %>% dplyr::filter(is.driver) %>% 
      tidyr::pivot_wider(names_from="sampleID", values_from="CCF")

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
    
    saveRDS(object=trees, file=paste0(out_dirname, "ctree.rds"))

    pdf(paste0(out_dirname, "plot_tree.pdf"))
    print(plot_tree)
    dev.off()
    """
}
