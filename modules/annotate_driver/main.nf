process ANNOTATE_DRIVER {
    publishDir params.publish_dir, mode: 'copy'

    input:

    tuple val(datasetID), val(patientID), val(sampleID), path(snv_RDS) 

    output:

    tuple val(datasetID), val(patientID), val(sampleID), path("$datasetID/$patientID/$sampleID/driver_annotation/*.rds"), emit: rds

    script:

    """
    #!/usr/bin/env Rscript

    library(tidyverse)

    res_dir = paste0("$datasetID", "/", "$patientID", "/", "$sampleID", "/")
    dir.create(res_dir, recursive = TRUE)

    SNV = readRDS("$snv_RDS")
    SNV = SNV[["$sampleID"]]
    SNV = SNV\$mutations

    drivers_table = read_tsv(file = '~/Downloads/2023-05-31_IntOGen-Drivers/Compendium_Cancer_Genes.tsv')

    x = SNV %>% 
      mutate(CANCER_TYPE = $cancer_type) %>% // we should definite it somewhere, maybe in the input sample sheet, maybe in place of dataset id
      dplyr::left_join(
        drivers_table %>% 
          dplyr::select(SYMBOL, CANCER_TYPE, CGC_CANCER_GENE) %>% 
          unique(),
        by = c('SYMBOL', 'CANCER_TYPE')) %>% 
      dplyr::mutate(is_driver = (CGC_CANCER_GENE & IMPACT %in% c('MODERATE', 'HIGH'))) %>% // this is our rule for driver assignment

    saveRDS(object = x, file = paste0(res_dir, "annotated_drivers.rds"))

    """
}
