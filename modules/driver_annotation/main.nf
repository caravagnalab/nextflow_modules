process ANNOTATE_DRIVER {
    publishDir params.publish_dir, mode: 'copy'

    input:

    tuple val(datasetID), val(patientID), val(sampleID), path(rds_list, stageAs: '*.rds') 

    output:

    tuple val(datasetID), val(patientID), val(sampleID), path("$datasetID/$patientID/$sampleID/driver_annotation/*.rds"), emit: rds

    script:

    """
    #!/usr/bin/env Rscript

    library(tidyverse)

    res_dir = paste0("$datasetID", "/", "$patientID", "/", "$sampleID", "/")
    dir.create(res_dir, recursive = TRUE)
    """
}