process MOBSTERh {
  publishDir (
    params.publish_dir,
    mode: "copy"
  )

  input:
    tuple val(patientID), val(sampleID), path(joint_table)

  output:
    //tuple val(patientID), path("$patientID/ctree/ctree_input_mobsterh.csv"), emit: ctree_input  // do not save or save inside mobster
    // tuple val(patientID), path("$patientID/*/mobster/*.rds")  // save also fits for each sample in mobster/sample_id/mobster.rds
    tuple val(patientID), val(sampleID), path("$patientID/$sampleID/mobster/*.rds"), emit: mobster_rds
    tuple val(patientID), val(sampleID), path("$outDir/mobster_joint*"), emit: mobster_joint
    //tuple val(patientID), path("$patientID/mobster/*.pdf"), emit: mobster_pdf

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
    outDir = "$patientID/$sampleID/mobster"
    new_joint = "$patientID/$sampleID/mobster/mobster_joint_table.tsv"
 

    """
    #!/usr/bin/env Rscript

    # Sys.setenv("VROOM_CONNECTION_SIZE"=99999999)

    library(mobster)
    library(dplyr)
    source("$moduleDir/getters.R")

    patientID = "$patientID"
    description = patientID

    ## Extract col names from joint
    out_dirname_mobsterh = paste0("$patientID","/","$sampleID","/mobster/")
    dir.create(out_dirname_mobsterh, recursive = TRUE)
    samples <-strsplit(x = "$sampleID", " ")%>% unlist()
    print(samples)
    
    orginal <- readRDS("$joint_table") %>% get_sample(sample = samples,which_obj = "orginal")
    joint_table = lapply(names(original), function(sample_name) CNAqc::Mutations(x=original[[sample_name]]) %>% dplyr::mutate(sample_id=sample_name)) %>% 
      dplyr::bind_rows()
    #read.csv("$joint_table", sep="\t") %>% filter(sample_id%in%samples) %>% 
    #  write.table(append = F,file = paste0(out_dirname_mobsterh,"joint_table.tsv"), quote = F,sep = "\t",row.names = F)

    
    # input_tab = read.csv("$outDir/joint_table.tsv", sep="\t") %>%
    input_tab = joint_table %>%
      #dplyr::rename(variantID = gene) %>%
      #dplyr::rename(is.driver = is_driver) %>%
      #dplyr::rename(tumour_content = purity) %>%
      dplyr::mutate(VAF=replace(VAF, VAF==0, 1e-7))
      #dplyr::rename(is_driver=is.driver) 
      #dplyr::rename(is_driver=is.driver, driver_label=variantID)

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
    fit = run_mobster_fit(inp_tb = input_tab, descr = description)
    best_fit = fit[["best"]]
    plot_fit = plot(best_fit)
    annotated_tab = best_fit[["data"]]
    sample_id = "$sampleID"
    write.table(x = annotated_tab,file = paste0("$outDir","/mobster_joint_",sample_id,".tsv"),append = F,quote = F,sep = "\t",row.names = F)
    saveRDS(object=fit, file=paste0(out_dirname_mobsterh, "mobsterh_fit.rds"))
    """
}
