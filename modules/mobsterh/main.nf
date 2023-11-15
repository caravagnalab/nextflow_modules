process MOBSTERh {
  publishDir params.publish_dir

  input:
    tuple val(patientID), val(timepointID), val(sampleID), path(joint_table)

  output:
    path("$patientID/$timepointID/$sampleID/*.rds")

  script:
    def args = task.ext.args ?: ''
    def K = args!='' && args.K ? "$args.K" : "1:3"
    def samples = args!='' && args.samples ? "$args.samples" : "5"
    def init = args!='' && args.init ? "$args.init" : "peaks"
    def tail = args!='' && args.tail ? "$args.tail" : "c(TRUE,FALSE)"
    def epsilon = args!='' && args.epsilon ? "$args.epsilon" : "1e-10"
    def maxIter = args!='' && args.maxIter ? "$args.maxIter" : "250"
    def fit_type = args!='' && args.fit_type ? "$args.fit_type" : "MM"
    def seed = args!='' && args.seed ? "$args.seed" : "12345"
    def model_selection = args!='' && args.model_selection ? "$args.model_selection" : "reICL"
    def trace = args!='' && args.trace ? "$args.trace" : "FALSE"
    def parallel = args!='' && args.parallel ? "$args.parallel" : "TRUE"
    def pi_cutoff = args!='' && args.pi_cutoff ? "$args.pi_cutoff" : "0.02"
    def n_cutoff = args!='' && args.n_cutoff ? "$args.n_cutoff" : "10"
    def auto_setup = args!='' && args.auto_setup ? "$args.auto_setup" : "NULL"
    def silent = args!='' && args.silent ? "$args.silent" : "FALSE"
    
    """
    #!/usr/bin/env Rscript

    # Sys.setenv("VROOM_CONNECTION_SIZE"=99999999)

    library(mobster)
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
    p = plot(best_fit) 

    dir.create(paste0("$patientID","/","$timepointID","/","$sampleID"), recursive = TRUE)
    saveRDS(object=fit, file=paste0("$patientID","/","$timepointID","/","$sampleID","/mobsterh.rds"))
    ggsave(filename = paste0("$patientID","/","$timepointID","/","$sampleID","/mobsterh.pdf"), plot=p)
    """
}
