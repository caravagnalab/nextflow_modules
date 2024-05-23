//
// Mutations extraction from mCNAqc
//

process RDS_PROCESSING {
    publishDir params.publish_dir, mode: 'copy'

    input:
      tuple val(datasetID), val(patientID),  path(join_cnaqc)

    output:
      tuple val(datasetID), val(patientID), path("$datasetID/$patientID/CNAqc2tsv/joint_table.tsv")
    
    script:

      def args                              = task.ext.args                                 ?: ''


    """
    #!/usr/bin/env Rscript

    library(dplyr)
    library(CNAqc)
    library(tidyr)

    source("$moduleDir/utils.R")

 
    multi_cnaqc = readRDS(file = "$join_cnaqc")
    mutations_multisample <- get_sample(m_cnaqc_obj = multi_cnaqc,sample = get_sample_name(multi_cnaqc),
                                    which_obj = "original")
    multisample_jointTable <- list()

    for (s in get_sample_name(multi_cnaqc)){
      purity <- mutations_multisample[[s]][["purity"]]
      multisample_joint[[s]] <- mutations_multisample[[s]][["mutations"]] %>%
        dplyr::mutate(purity=purity)
      }
    
    joint_table <- bind_rows(multisample_jointTable)

    write.table(joint_table, file = "$datasetID/$patientID/CNAqc2tsv/joint_table.tsv", append = F, quote = F, sep = "\t", row.names = FALSE)

    """
}
