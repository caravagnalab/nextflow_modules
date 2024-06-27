//
// Mutations extraction from mCNAqc
//

process RDS_PROCESSING {
    publishDir params.publish_dir, mode: 'copy'

    input:
      tuple val(datasetID), val(patientID),  val(sampleID), path(join_cnaqc)

    output:
      tuple val(datasetID), val(patientID), val(sampleID), path("formatter/CNAqc2tsv/$datasetID/$patientID/joint_table_mut.tsv"), emit: tsv
      tuple val(datasetID), val(patientID), val(sampleID), path("formatter/CNAqc2tsv/$datasetID/$patientID/joint_table_cnv.tsv"), emit: tsv
    
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
    multisample <- get_sample(m_cnaqc_obj = multi_cnaqc,sample = get_sample_name(multi_cnaqc),
                                    which_obj = "original")


    # Extract mutations data

    multisampleMut_jointTable <- list()
    

    for (s in get_sample_name(multi_cnaqc)){
      purity <- mutations_multisample[[s]][["purity"]]
      multisampleMut_jointTable[[s]] <- multisample[[s]][["mutations"]] %>%
        dplyr::mutate(purity=purity)
      }
    
    joint_table_mut <- bind_rows(multisampleMut_jointTable)

    write.table(joint_table_mut, file = paste0(res_dir,"joint_table_mut.tsv"), append = F, quote = F, sep = "\t", row.names = FALSE)

    # Extract CNA data
    
    multisampleCNA_jointTable <- list()

    for (s in get_sample_name(multi_cnaqc)){
      multisampleCNA_jointTable[[s]] <- multisample[[s]][["cna"]] %>% 
      dplyr::mutate(sample_id = s)
    }
    
    joint_table_cna <- bind_rows(multisampleCNA_jointTable)

    write.table(joint_table_cna, file = paste0(res_dir,"joint_table_cna.tsv"), append = F, quote = F, sep = "\t", row.names = FALSE)
    
    """
}
