process VCF_PROCESSING {
    publishDir params.publish_dir, mode: 'copy'

    input:
     tuple val(datasetID), val(patientID), val(sampleID), path(vcfFile), path(vcfTbi)

    output:
     tuple val(datasetID), val(patientID), val(sampleID), path("$datasetID/$patientID/$sampleID/vcf2CNAqc/*.rds"), emit: rds 

    script:

    """
    #!/usr/bin/env Rscript 
    
    library(tidyverse)
    library(vcfR)
    
    res_dir = paste0("$datasetID", "/", "$patientID", "/", "$sampleID", "/vcf2CNAqc/")
    dir.create(res_dir, recursive = TRUE)

    source(paste0("$moduleDir", '/parser_vcf.R'))

    # Read vcf file
    vcf = vcfR::read.vcfR("$vcfFile")

    # Check from which caller the .vcf has been produced
    source = queryMETA(vcf, element = 'source')[[1]]

    if (grepl(pattern = 'Mutect', x = source)){
        calls = parse_Mutect(vcf, sample_id = "$sampleID")
        
    } else if {
        calls = NA
    }

    saveRDS(object = calls, file = paste0(res_dir, "VCF.rds"))
    """
}
