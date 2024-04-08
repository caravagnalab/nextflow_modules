process GET_POSITIONS {
    publishDir params.publish_dir, mode: 'copy'

    input:

    tuple val(datasetID), val(patientID), val(sampleID), path(rds_list, stageAs: '*.rds') 

    output:

    tuple val(datasetID), val(patientID), val(sampleID), path("$datasetID/$patientID/*/mpileup/*.bed"), emit: bed
    tuple val(datasetID), val(patientID), val(sampleID), path("$datasetID/$patientID/mpileup/*.bed"), emit: pos

    script:

    """
    #!/usr/bin/env Rscript

    library(tidyverse)

    res_dir = paste0("$datasetID", "/", "$patientID", "/")
    dir.create(res_dir, recursive = TRUE)
    
    samples = substr("$sampleID", 2, nchar("$sampleID")-1)
    samples = strsplit(samples, ", ")[[1]]
    #samples = strsplit(x = "$sampleID", " ")%>% unlist()

    positions = lapply(strsplit("$rds_list", " ")[[1]], FUN = function(rds){
        df = readRDS(rds)
        df = df[[1]]\$mutations %>% dplyr::mutate(id = paste(chr, from, to, sep = ":")) %>% dplyr::select(chr, from, to, ref, alt, id)
    }) 

    all_positions = positions %>% dplyr::bind_rows() %>% dplyr::pull(id) %>% unlist() %>% unique()
    all = positions %>% dplyr::bind_rows() %>% dplyr::distinct() %>% dplyr::select(-id)

    dir.create(paste0(res_dir, "mpileup/"), recursive = TRUE)
    write.table(file = paste0(res_dir, 'mpileup/all_positions.bed'), all, quote = F, sep = "\t", row.names = F, col.names = T)

    missed = lapply(seq(1,length(positions)), FUN = function(s){
        all_positions[!(all_positions %in% positions[[s]]\$id)]
    }) 

    for (i in seq(1,length(missed))){
        sample = samples[[i]]
        df = dplyr::tibble(id = missed[[i]]) %>% 
                tidyr::separate(id,into = c('chr', 'from', 'to'), sep = ':') %>% 
                dplyr::filter(chr %in% c(paste0('chr', seq(1,22)), 'chrX', 'chrY'))

	dir.create(paste0(res_dir, sample, "/mpileup/"), recursive = TRUE)                
        write.table(file = paste0(res_dir, sample, "/mpileup/positions_missing", ".bed"), df, quote = F, sep = "\t", row.names = F, col.names = F)
    }
    """
}
