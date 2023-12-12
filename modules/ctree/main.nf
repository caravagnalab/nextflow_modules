process CTREE {
  publishDir params.publish_dir, mode: 'copy'

  input:
    tuple val(patientID), val(timepointID), val(sampleID), path(ccf_table)

  output:
    tuple path("$patientID/$timepointID/$sampleID/*.rds"), path("$patientID/$timepointID/$sampleID/*.pdf")

  script:
    def args = task.ext.args ?: ''
    def sspace_cutoff = args!='' && args.sspace_cutoff ? "$args.sspace_cutoff" : ""
    def n_sampling = args!='' && args.n_sampling ? "$args.n_sampling" : ""
    def store_max = args!='' && args.store_max ? "$args.store_max" : ""

    """
    #!/usr/bin/env Rscript

    library(ctree)
    library(dplyr)

    input_tab = read.csv("$ccf_table")

    # the CCF table must report CCF values for each cluster and sample
    # cluster | nMuts | is.driver | is.clonal | sample1 | sample2 | ...
    CCF_table = input_tab %>% 
        # remove mutations in the tail
        dplyr::filter(cluster != "Tail")

    # the driver table must contain patient and variant IDs and report clonality and driver status
    # patientID | variantID | is.driver | is.clonal | cluster | sample1 | sample2 | ...
    drivers_table = input_tab # ... reformat

    samples = c("$sampleID")  # if multisample, this is a list
    patient = "$patientID"


    ## from MOBSTER object
    cluster_table = mobster_fit[["Clusters"]] %>% 
        dplyr::filter(cluster != "Tail", type == "Mean") %>% 
        dplyr::select(cluster, fit.value) %>% 
        rename(R1 = fit.value)
    cluster_table[["nMuts"]] = mobster_fit[["N.k"]][cluster_table[["cluster"]]]
    cluster_table[["is.clonal"]] = FALSE
    cluster_table[["is.clonal"]][which.max(cluster_table[["R1"]])] = TRUE

    drivers_collapse = mobster_fit[["data"]] %>% 
        dplyr::filter(is_driver) %>% 
        dplyr::pull(cluster) %>% unique
    cluster_table[["is.driver"]] = FALSE
    cluster_table[["is.driver"]][which(cluster_table[["cluster"]] %in% drivers_collapse)] = TRUE
    drivers_table = mobster_fit[["data"]] %>% dplyr::as_tibble() %>% 
        dplyr::filter(is_driver) %>% 
        dplyr::rename(variantID = driver_label, is.driver = is_driver) %>% 
        dplyr::mutate(patientID = patientID, R1 = VAF)
    drivers_table[["is.clonal"]] = FALSE
    drivers_table[["is.clonal"]][which(drivers_table[["cluster"]] == cluster_table %>% 
        dplyr::filter(is.clonal) %>% dplyr::pull(cluster))] = TRUE
    drivers_table = drivers_table %>% dplyr::select(patientID, variantID, is.driver, is.clonal, cluster, R1, dplyr::everything())

    trees = ctrees(CCF_clusters = cluster_table,
                  drivers = drivers_table,
                  samples = samples,
                  patient = patient,
                  sspace.cutoff = as.integer("$sspace_cutoff"),
                  n.sampling = as.integer("$n_sampling"),
                  store.max = as.integer("$store_max"))

    # plot the best tree
    plot_tree = plot(trees[[1]]) 

    dir.create(paste0("$patientID","/","$timepointID","/","$sampleID"), recursive = TRUE)
    saveRDS(object=trees, file=paste0("$patientID","/","$timepointID","/","$sampleID","/ctree.rds"))
    pdf(paste0("$patientID","/","$timepointID","/","$sampleID","/plot_tree.pdf"))
    print(plot_tree)
    dev.off()
    """
}
