process VCF_PROCESSING {
    publishDir params.publish_dir, mode: 'copy'

    input:
     tuple val(datasetID), val(patientID), val(sampleID), path(vcfFile)

    output:
     tuple val(datasetID), val(patientID), val(sampleID), path("$datasetID/$patientID/$sampleID/vcf2CNAqc/*.rds"), emit: rds 

    script:
        def args              = task.ext.args                         ?: ''
        def filter_mutations  = args!='' && args.filter_mutations     ?  "$args.filter_mutations" : "FALSE"

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
    source = vcfR::queryMETA(vcf, element = 'source')[[1]]

    if (grepl(pattern = 'Mutect', x = source)){
        calls = parse_Mutect(vcf, sample_id = "$sampleID", filter_mutations = "$filter_mutations")
        
    } else if (grepl(pattern = 'Strelka', x = source)){
        calls = parse_Strelka(vcf, sample_id = "$sampleID", filter_mutations = "$filter_mutations")
    
    } else if (grepl(pattern = 'Platypus', x = source)){
        calls = parse_Platypus(vcf, sample_id = "$sampleID", filter_mutations = "$filter_mutations")

    } else if (grepl(pattern = 'freeBayes', x = source)){
        calls = parse_Freebayes(vcf, sample_id = "$sampleID", filter_mutations = "$filter_mutations") 

    } else {
        stop('Variant Caller not supported.')
    }

    saveRDS(object = calls, file = paste0(res_dir, "VCF.rds"))
    """
}
