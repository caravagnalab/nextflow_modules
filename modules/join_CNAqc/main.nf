process JOIN_CNAQC {
  publishDir params.publish_dir, mode: 'copy'

  input:
    
    tuple val(datasetID), val(patientID), val(sampleID), path(rds_list, stageAs: '*.rds')
  
  output:

    tuple val(datasetID), val(patientID), path("$datasetID/$patientID/join_CNAqc/*.rds"), emit: rds

  script:
    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    library(CNAqc)
    
    res_dir = paste0("$datasetID", "/", "$patientID", "/join_CNAqc/")
    dir.create(res_dir, recursive = TRUE)
    
    samples = substr("$sampleID", 2, nchar("$sampleID")-1)
    samples = strsplit(samples, ", ")[[1]]

    source(paste0("$moduleDir", '/join_CNAqc.R'))

    result = lapply(strsplit("$rds_list", " ")[[1]], FUN = function(file){
             readRDS(file)
             }) 
    names(result) = samples
    
    for (name in names(result)){
      result[[name]]\$mutations = result[[name]]\$mutations %>% dplyr::rename(Indiv = sample)
    }
    
    out = multisample_init(result)
    saveRDS(object = out, file = paste0(res_dir, "prova_multi_cnaqc.rds"))
    """
}
