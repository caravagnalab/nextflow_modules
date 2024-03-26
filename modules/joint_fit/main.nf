process JOINT_FIT {
  publishDir (
    params.publish_dir,
    mode: "copy"
  )

  input:
    tuple val(patientID), val(sampleID), path(mobster_joint)

  output:
    //tuple val(patientID), path("$patientID/ctree/ctree_input_mobsterh.csv"), emit: ctree_input  // do not save or save inside mobster
    // tuple val(patientID), path("$patientID/*/mobster/*.rds")  // save also fits for each sample in mobster/sample_id/mobster.rds
    tuple val(patientID), val(sampleID), path("$patientID/mobster/annotated_joint_table.tsv"), emit: annotated_joint
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
    def tools = args!="" && args.tools ? "$args.tools" : ""

    """
    #!/usr/bin/env Rscript

    # Sys.setenv("VROOM_CONNECTION_SIZE"=99999999)
    library(plyr)
    library(dplyr)
    patientID = "$patientID"
    description = patientID
    dir.create("$patientID/mobster/", recursive = TRUE) 
    tsv_files = strsplit("$mobster_joint", " ")[[1]]
    tsvs <- list()
    for (i in seq_along(tsv_files)){
       tsvs[[i]]<-read.table(file = tsv_files[i],header = T,sep = "\t")
    }
    t = plyr::ldply(tsvs, rbind)
    t["cluster"] %>% 
     apply(MARGIN=1, FUN = function(w) {any(w=="Tail", na.rm=TRUE)}) %>% table()
    non_tail = t %>% filter(cluster!="Tail") %>% rename(mobster_cluster_id = cluster)
    write.table(x = non_tail,file = "$patientID/mobster/annotated_joint_table.tsv",append = F,quote = F,sep = "\t",row.names = F)
    
    """
}
