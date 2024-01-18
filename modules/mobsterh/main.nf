process MOBSTERh {
  publishDir params.publish_dir, mode: "copy"

  input:
    tuple val(patientID), path(joint_table)

  output:
    tuple val(patientID), path("$patientID/ctree/ctree_input_mobsterh.csv"), emit: ctree_input  // do not save or save inside mobster
    // tuple val(patientID), path("$patientID/*/mobster/*.rds")  // save also fits for each sample in mobster/sample_id/mobster.rds
    tuple val(patientID), path("$patientID/mobster/*.rds"), emit: mobster_rds  // save also fits for each sample in mobster/sample_id/mobster.rds
    tuple val(patientID), path("$patientID/mobster/*.pdf"), emit: mobster_pdf

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

    patientID = "$patientID"
    description = patientID
    input_tab = read.csv("$joint_table", sep="\t") %>%
      dplyr::filter(patientID==patientID) %>% 
      dplyr::mutate(DP=ref_counts+alt_counts, 
                    VAF=alt_counts/DP) %>%
      # dplyr::filter(VAF <= 1) %>% 
      dplyr::mutate(VAF=replace(VAF, VAF==0, 1e-7)) %>% 
      dplyr::rename(is_driver=is.driver, driver_label=variantID)

    multisample = length(unique(input_tab[["sample_id"]])) > 1

    print(multisample)

    run_mobster_fit = function(inp_tb, descr) {
      mobster_fit(x = inp_tb,
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
                  description = descr)
    }

    ## If single sample
    if (!multisample) {
      fit = run_mobster_fit(inp_tb = input_tab, descr = description)

      best_fit = fit[["best"]]
      plot_fit = plot(best_fit)

      cluster_table = best_fit[["Clusters"]] %>% 
        dplyr::filter(cluster != "Tail", type == "Mean") %>%
        dplyr::select(cluster, fit.value) %>% 
        dplyr::rename(R1 = fit.value) %>% 
        dplyr::mutate(patientID = patientID)
      cluster_table[["nMuts"]] = best_fit[["N.k"]][cluster_table[["cluster"]]]
      cluster_table[["is.clonal"]] = FALSE
      cluster_table[["is.clonal"]][which.max(cluster_table[["R1"]])] = TRUE

      drivers_collapse = best_fit[["data"]] %>% dplyr::filter(is_driver) %>% 
        dplyr::pull(cluster) %>% unique

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

      sample_names = "R1"

    } else {
      ## If multisample
      library(VIBER)

      samples = unique(input_tab[["sample_id"]])
      fits = lapply(samples,
                    function(s) {
                      input_tab %>%
                        dplyr::filter(sample_id == s) %>%
                        # dplyr::filter(VAF > 0.05) %>%
                        run_mobster_fit(descr = paste0(description, "_", s))
                    }) %>% setNames(samples)
      
      myjoin = function(x, y) dplyr::full_join(x, y, by="mutation_id")

      cols_to_remove = c("patientID", "sample_id", "driver_label", "is_driver")
      assignments = lapply(names(fits), function(fi) {
        D = fits[[fi]][["best"]][["data"]] %>% 
          dplyr::select(-dplyr::all_of(cols_to_remove)) %>% 
          dplyr::select(mutation_id, dplyr::everything())
        colnames(D)[2:ncol(D)] = paste0(colnames(D)[2:ncol(D)], ".", fi)
        D
      }) %>% Reduce(myjoin, .)
      
      # Find tail mutations
      ## check which mutations are removed -> leave as choice if remove muts in tail for ALL or at least one sample
      ## this removes mutations that are in the tail in at least one sample
      assignments[["tail"]] = assignments %>%
        dplyr::select(dplyr::starts_with("cluster")) %>%
        apply(MARGIN=1, FUN = function(w) {any(w=="Tail", na.rm=TRUE)})
      
      # Non-tail mutations
      non_tail = assignments %>% dplyr::filter(!tail) %>% dplyr::pull(mutation_id)
      
      # # Read counts
      # reads_data = lapply(samples, function(s) {
      #   input_tab %>% dplyr::filter(sample_id == s)
      # }) %>% Reduce(myjoin, .)

      reads_data = input_tab %>% 
        tidyr::pivot_wider(names_from="sample_id",
                          values_from=c("ref_counts","alt_counts","normal_cn",
                                        "major_cn","minor_cn","tumour_content",
                                        "DP","VAF"), names_sep=".")
      
      dp = reads_data %>% dplyr::filter(mutation_id %in% non_tail) %>% 
        dplyr::select(dplyr::starts_with("DP")) %>% 
        dplyr::mutate(dplyr::across(.cols=dplyr::everything(), 
                                    .fns=function(x) replace(x, is.na(x), 0))) %>% 
        dplyr::rename_all(function(x) stringr::str_remove_all(x,"DP."))
      nv = reads_data %>% dplyr::filter(mutation_id %in% non_tail) %>% 
        dplyr::select(dplyr::starts_with("alt_counts")) %>% 
        dplyr::mutate(dplyr::across(.cols=dplyr::everything(), 
                                    .fns=function(x) replace(x, is.na(x), 0))) %>% 
        dplyr::rename_all(function(x) stringr::str_remove_all(x,"alt_counts."))
      
      # Fit and plot
      viber_K = eval(parse(text="$K")) 
      viber_K[which.min(viber_K)] = 2
      fit = VIBER::variational_fit(nv, dp, 
                                   K=unique(viber_K), 
                                   data=reads_data %>% dplyr::filter(mutation_id %in% non_tail))

      best_fit = fit
      plot_fit = plot(best_fit)

      ## generate cluster table for ctree
      pi = best_fit[["pi_k"]][((best_fit[["N"]] * best_fit[["pi_k"]]) %>% round) > 0]
      theta = best_fit[["theta_k"]][ , names(pi), drop = T]
      
      # Get clusters table: cluster and fit
      cluster_table = data.frame(cluster=colnames(theta), stringsAsFactors=FALSE) %>% tibble::as_tibble()
      cluster_table = bind_cols(cluster_table, t(theta) %>% as_tibble)
      
      cluster_table[["nMuts"]] = table(best_fit[["labels"]])[cluster_table[["cluster"]]] %>% as.vector()
      
      # Clonality status - maximum fit is the clonal
      clonal_cluster = apply(theta, 1, which.max)
      clonal_cluster = colnames(theta)[clonal_cluster]
      clonal_cluster = which.max(table(clonal_cluster)) %>% names
      
      cluster_table[["is.clonal"]] = FALSE
      cluster_table[["is.clonal"]][cluster_table[["cluster"]] %in% clonal_cluster] = TRUE
      cluster_table = cluster_table %>% dplyr::mutate(patientID=patientID)
      
      best_fit[["data"]][["cluster"]] = paste(unlist(best_fit[["labels"]]))
      
      drivers_collapse = best_fit[["data"]] %>%
        dplyr::filter(is_driver) %>%
        dplyr::pull(cluster) %>%
        unique()
      
      cluster_table[["is.driver"]] = FALSE
      cluster_table[["is.driver"]][which(cluster_table[["cluster"]] %in% drivers_collapse)] = TRUE
      
      # Create drivers table
      cx = best_fit[["x"]] %>% dplyr::select(-dplyr::starts_with("cluster"))
      cy = best_fit[["y"]] %>% dplyr::select(-dplyr::starts_with("cluster"))
      vaf_table = cx/cy
      
      drivers_table = best_fit[["data"]] %>%
        dplyr::filter(is_driver) %>% 
        dplyr::rename(variantID=driver_label, is.driver=is_driver)
      # drivers_table = dplyr::bind_cols(drivers_table, vaf_table[which(best_fit[["data"]][["is.driver"]]), , drop = F])
      
      drivers_table[["is.clonal"]] = FALSE
      drivers_table[["is.clonal"]][which(
        drivers_table[["cluster"]] == cluster_table %>% dplyr::filter(is.clonal) %>% dplyr::pull(cluster)
      )] = TRUE
      
      drivers_table[["patientID"]] = patientID

      sample_names = colnames(cx)

      drivers_table = drivers_table %>%
        dplyr::select(patientID, variantID, is.driver, is.clonal, cluster, dplyr::everything())

    }

    # generate output for ctree
    ctree_input = dplyr::full_join(drivers_table, cluster_table, 
                    by=c("patientID", "is.driver", "is.clonal", "cluster")) %>% 
      tidyr::pivot_longer(cols=sample_names, names_to="sampleID", values_to="CCF")

    ## save rds and plots
    out_dirname_mobsterh = paste0("$patientID", "/mobster/")
    out_dirname_ctree = paste0("$patientID","/ctree/")

    dir.create(out_dirname_mobsterh, recursive = TRUE)
    dir.create(out_dirname_ctree, recursive = TRUE)

    saveRDS(object=fit, file=paste0(out_dirname_mobsterh, "mobsterh_fit.rds"))
    saveRDS(object=plot_fit, file=paste0(out_dirname_mobsterh, "mobsterh_plot.rds"))
    
    pdf(paste0(out_dirname_mobsterh, "mobsterh.pdf"))
    print(plot_fit)
    dev.off()

    write.csv(x=ctree_input, file=paste0(out_dirname_ctree, "ctree_input_mobsterh.csv"), row.names=F)
    """
}
