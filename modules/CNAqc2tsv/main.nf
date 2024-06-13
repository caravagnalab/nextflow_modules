//
// Mutations extraction from mCNAqc
//

process RDS_PROCESSING {
    publishDir params.publish_dir, mode: 'copy'

    input:
      tuple val(datasetID), val(patientID),  val(sampleID), path(join_cnaqc)

    output:
      tuple val(datasetID), val(patientID), val(sampleID), path("formatter/CNAqc2tsv/$datasetID/$patientID/joint_table.tsv"), emit: tsv
    
    script:

      def args                              = task.ext.args                                 ?: ''


    """
    #!/usr/bin/env Rscript

    library(dplyr)
    library(CNAqc)
    library(tidyr)

    source("$moduleDir/utils.R")
    
    res_dir = paste0("formatter/CNAqc2tsv/","$datasetID", "/", "$patientID/")
    dir.create(res_dir, recursive = TRUE)
 
    multi_cnaqc = readRDS(file = "$join_cnaqc")
    mutations_multisample <- get_sample(m_cnaqc_obj = multi_cnaqc,sample = get_sample_name(multi_cnaqc),
                                    which_obj = "original")
    multisample_jointTable <- list()

    for (s in get_sample_name(multi_cnaqc)){
      purity <- mutations_multisample[[s]][["purity"]]
      multisample_jointTable[[s]] <- mutations_multisample[[s]][["mutations"]] %>%
        dplyr::mutate(purity=purity)
      }
    
    joint_table <- bind_rows(multisample_jointTable)

    write.table(joint_table, file = paste0(res_dir,"joint_table.tsv"), append = F, quote = F, sep = "\t", row.names = FALSE)
    
    """
}
