process CNA_PROCESSING {
    publishDir params.publish_dir, mode: 'copy'

    input:
     tuple val(datasetID), val(patientID), val(sampleID), path(cnaDir)

    output:
     tuple val(datasetID), val(patientID), val(sampleID), path("$datasetID/$patientID/$sampleID/cna2CNAqc/CNA.rds"), emit: rds

    script:
    
      def args              = task.ext.args                         ?: ''
      def cna_caller        = args!='' && args.cna_caller           ? "$args.cna_caller" : "sequenza"

    """
    #!/usr/bin/env Rscript 
    
    library(tidyverse)
    library(vcfR)

    res_dir = paste0("$datasetID", "/", "$patientID", "/", "$sampleID", "/cna2CNAqc/")
    dir.create(res_dir, recursive = TRUE)

    if ("$cna_caller" == "sequenza"){
        load_SQ_output = function(sample, run){
              
              segments_file = paste0(run, '/', sample, '_segments.txt')
              purity_file = paste0(run, '/', sample, '_confints_CP.txt')

              segments = readr::read_tsv(segments_file, col_types = readr::cols()) %>%
                dplyr::rename(
                  chr = chromosome,
                  from = start.pos,
                  to = end.pos,
                  Major = A,
                  minor = B
                ) %>%
                dplyr::select(chr, from, to, Major, minor, dplyr::everything())

              solutions = readr::read_tsv(purity_file, col_types = readr::cols())

              purity = solutions\$cellularity[2]
              ploidy = solutions\$ploidy.estimate[2]
              
              return(list(
                  segments = segments,
                  purity = purity,
                  ploidy = ploidy
                )
              )
            }

        CNA = load_SQ_output("$patientID", "$cnaDir")
    
    } else { stop("Not valid CNA caller!") }

    saveRDS(object = CNA, file = paste0(res_dir, "CNA.rds"))
    """
}