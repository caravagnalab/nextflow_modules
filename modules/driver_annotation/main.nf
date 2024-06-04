process ANNOTATE_DRIVER {
    publishDir params.publish_dir, mode: 'copy'

    input:

    tuple val(datasetID), val(patientID), val(sampleID), path(snv_RDS)
    val(cancer_type) 

    output:

    tuple val(datasetID), val(patientID), val(sampleID), path("DriverAnnotation/$datasetID/$patientID/$sampleID/*.rds"), emit: rds

    script:

    """
    #!/usr/bin/env Rscript

    library(dplyr)
    library(readr)

    res_dir = paste0("DriverAnnotation/", "$datasetID", "/", "$patientID", "/", "$sampleID", "/")
    dir.create(res_dir, recursive = TRUE)

    data = readRDS("$snv_RDS")
    SNV = data[["$sampleID"]]
    SNV = SNV\$mutations

    drivers_table = readr::read_tsv(file = "$params.drivers_table") 
    
    if("$cancer_type" == 'PANCANCER'){
      drivers_table = drivers_table %>% 
        dplyr::group_by(SYMBOL) %>% 
        dplyr::reframe(CGC_CANCER_GENE = any(CGC_CANCER_GENE), dplyr::across(dplyr::everything())) %>% 
        dplyr::filter(CGC_CANCER_GENE) %>% 
        dplyr::mutate(CANCER_TYPE = 'PANCANCER')
    } 

    drivers_table = drivers_table %>% 
        dplyr::select(SYMBOL, CANCER_TYPE, CGC_CANCER_GENE) %>% 
        unique()


    x = SNV %>% 
      dplyr::mutate(CANCER_TYPE = "$cancer_type") %>%
      dplyr::left_join(
        drivers_table,
        by = c('SYMBOL', 'CANCER_TYPE')
      ) %>% 
      dplyr::mutate(
          is_driver = (CGC_CANCER_GENE & IMPACT %in% c('MODERATE', 'HIGH')),
          driverl_label = paste(SYMBOL, HGVSp)
      )

    new_data = list()
    new_data[["$sampleID"]]\$mutations = x
    new_data[["$sampleID"]]\$sample = data[["$sampleID"]]\$sample
    saveRDS(object = new_data, file = paste0(res_dir, "annotated_drivers.rds"))

    """
}