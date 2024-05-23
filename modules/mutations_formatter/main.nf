//
// Mutations extraction from mCNAqc
//

process MUTATIONS_FORMATTER {
    publishDir params.publish_dir, mode: 'copy'

    input:
      tuple val(datasetID), path(join_cnaqc)

    output:
      tuple val(datasetID), path("$datasetID/mutations_extraction/multisample_mutations.tsv")
    
    script:

      def args                              = task.ext.args                                 ?: ''


    """
    #!/usr/bin/env Rscript

    library("tidyverse")
    library("CNAqc")

    mutations_getter = paste0("res_mutGetter/")
    dir.create(mutations_getter, recursive = TRUE)

    #Input dataset : mCNAqc object
    multi_cnaqc <- readRDS('$join_cnaqc')

    #Extract mutations data
    mutations_multisample = CNAqc::get_sample(m_cnaqc_obj = multi_cnaqc,
                 sample = CNAqc::get_sample_name(multi_cnaqc),
                 which_obj = "original")

    mutations_multisample <- lapply(mutations_mcnaqc, function(x) {CNAqc::Mutations(x)}) %>% dplyr::bind_rows()

    write.table(mutations_multisample, file = "$datasetID/mutations_getter/mutations_multisample.tsv", row.names=TRUE, sep="\t")

    """
}
