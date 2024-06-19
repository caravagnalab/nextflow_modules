process JOIN_CNAQC {
  publishDir params.publish_dir, mode: 'copy'

  input:
    
    tuple val(datasetID), val(patientID), val(sampleID), path(rds_list, stageAs: '*.rds')
  
  output:

    tuple val(datasetID), val(patientID), val (sampleID), path("QC/join_CNAqc/$datasetID/$patientID/*.rds"), emit: rds

  script:

    def args                                    = task.ext.args                                             ?: ''
    def qc_filter                               = args!='' && args.qc_filter                                ?  "$args.qc_filter" : ""
    def keep_original                           = args!='' && args.keep_original                            ?  "$args.keep_original" : ""

    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    library(CNAqc)
    
    res_dir = paste0("QC/join_CNAqc/", "$datasetID", "/", "$patientID", "/")
    dir.create(res_dir, recursive = TRUE)
    
    samples = substr("$sampleID", 2, nchar("$sampleID")-1)
    samples = strsplit(samples, ", ")[[1]]

    result = lapply(strsplit("$rds_list", " ")[[1]], FUN = function(file){
             readRDS(file)
             }) 
    names(result) = samples
    
    for (name in names(result)){
      result[[name]]\$mutations = result[[name]]\$mutations %>% dplyr::rename(Indiv = sample)
    }
    
    out = CNAqc::multisample_init(result, 
                            QC_filter = as.logical("$qc_filter"), 
                            keep_original = as.logical("$keep_original"), 
                            discard_private = FALSE)
    saveRDS(object = out, file = paste0(res_dir, "multi_cnaqc.rds"))
    """
}
