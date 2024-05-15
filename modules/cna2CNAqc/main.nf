process CNA_PROCESSING {
    publishDir params.publish_dir, mode: 'copy'

    input:
     tuple val(datasetID), val(patientID), val(sampleID), path(cnaFile), val(caller)

    output:
     tuple val(datasetID), val(patientID), val(sampleID), path("$datasetID/$patientID/$sampleID/cna2CNAqc/CNA.rds"), emit: rds

    script:
    
    """
    #!/usr/bin/env Rscript 
  
    res_dir = paste0("$datasetID", "/", "$patientID", "/", "$sampleID", "/cna2CNAqc/")
    dir.create(res_dir, recursive = TRUE)

    source(paste0("$moduleDir", '/parser_CNA.R'))

    if ("$caller" == 'sequenza'){
      CNA = parse_Sequenza("$patientID", "$cnaFile")
    } else if {
      pass
    }

    saveRDS(object = CNA, file = paste0(res_dir, "CNA.rds"))
    """
}
