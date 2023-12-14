process MOBSTERh {
  publishDir params.publish_dir, mode: 'copy'

  input:
    tuple val(patientID), val(timepointID), val(sampleID), path(joint_table)

  output:
    tuple val(patientID), 
          path("$patientID/ctree/ctree_input.csv"),
          path("$patientID/$timepointID/$sampleID/*.rds"),
          path("$patientID/$timepointID/$sampleID/*.pdf")

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

    """
    #!/usr/bin/env Rscript

    # Sys.setenv("VROOM_CONNECTION_SIZE"=99999999)

    library(mobster)
    library(dplyr)

    description = paste("$patientID", "$timepointID", "$sampleID", sep="_")
    input_tab = read.csv("$joint_table")

    fit = mobster_fit(x = input_tab,
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
                      description = description)

    best_fit = fit[["best"]]
    plot_fit = plot(best_fit) 

    # generate output for ctree
    patientID = "$patientID"
    cluster_table = best_fit[["Clusters"]] %>% 
      dplyr::filter(cluster != "Tail", type == "Mean") %>%
      dplyr::select(cluster, fit.value) %>% 
      dplyr::rename(R1 = fit.value) %>% 
      dplyr::mutate(patientID = patientID)
    cluster_table[["nMuts"]] = best_fit[["N.k"]][cluster_table[["cluster"]]]
    cluster_table[["is.clonal"]] = FALSE
    cluster_table[["is.clonal"]][which.max(cluster_table[["R1"]])] = TRUE

    drivers_collapse = best_fit[["data"]] %>% dplyr::filter(is_driver) %>% 
      pull(cluster) %>% unique

    cluster_table[["is.driver"]] = FALSE
    cluster_table[["is.driver"]][which(cluster_table[["cluster"]] %in% drivers_collapse)] = TRUE
    drivers_table = best_fit[["data"]] %>% as_tibble() %>% dplyr::filter(is_driver) %>% 
      dplyr::rename(variantID = driver_label, is.driver = is_driver) %>% 
      dplyr::mutate(patientID = patientID, R1 = VAF) %>% 
      dplyr::select(-R1, -VAF)
    drivers_table[["is.clonal"]] = FALSE
    drivers_table[["is.clonal"]][which(drivers_table[["cluster"]] == cluster_table %>% 
                                    dplyr::filter(is.clonal) %>% dplyr::pull(cluster))] = TRUE
    drivers_table = drivers_table %>% 
      dplyr::select(patientID, variantID, is.driver, is.clonal, 
                    cluster, dplyr::everything())

    ctree_input = dplyr::full_join(drivers_table, cluster_table, 
                    by=c("patientID", "is.driver", "is.clonal", "cluster")) %>% 
      tidyr::pivot_longer(cols="R1", names_to="sampleID", values_to="CCF")

    # save rds and plots
    out_dirname_mobsterh = paste0("$patientID","/","$timepointID","/","$sampleID", "/")
    out_dirname_ctree = paste0("$patientID","/ctree/")

    dir.create(out_dirname_mobsterh, recursive = TRUE)
    dir.create(out_dirname_ctree, recursive = TRUE)

    saveRDS(object=fit, file=paste0(out_dirname_mobsterh, "mobsterh.rds"))
    
    pdf(paste0(out_dirname_mobsterh, "mobsterh.pdf"))
    print(plot_fit)
    dev.off()

    write.csv(x=ctree_input, file=paste0(out_dirname_ctree, "ctree_input.csv"), row.names=F)
    """
}
