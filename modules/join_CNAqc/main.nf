process JOIN_CNAQC {
  publishDir params.publish_dir, mode: 'copy'

  input:
    
    tuple val(datasetID), val(patientID), val(sampleID), path(rds_list, stageAs: '*.rds')
  
  output:

    tuple val(datasetID), val(patientID), val(sampleID), path("$datasetID/$patientID/join_CNAqc/*.rds"), emit: rds
    tuple val(datasetID), val(patientID), val(sampleID), path("$datasetID/$patientID/join_CNAqc/mut_join_table.tsv"), emit: mut_tsv
    tuple val(datasetID), val(patientID), val(sampleID), path("$datasetID/$patientID/join_CNAqc/cna_join_table.tsv"), emit: cna_tsv

  script:


    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    library(CNAqc)
    
    res_dir = paste0("$datasetID", "/", "$patientID", "/join_CNAqc/")
    dir.create(res_dir, recursive = TRUE)
    
    samples = substr("$sampleID", 2, nchar("$sampleID")-1)
    samples = strsplit(samples, ", ")[[1]]

    # Functions for JOIN_TABLE
    multisample_init = function(cnaqc_objs, cna_type = "clonal") {
        lapply(cnaqc_objs, function(x) {
            stopifnot(inherits(x, "cnaqc"))
        })
        
        multi_input = prepare_input_data_multiple(cnaqc_objs, cna_type)
        class(multi_input) <- "m_cnaqc"
        
        return(multi_input) 
    }  

    prepare_input_data_multiple = function(cnaqc_objs, cna_type) {
        #check if you are passing a list of cnaqc objs
        lapply(cnaqc_objs, function(x) {
            stopifnot(inherits(x, "cnaqc"))
        })
        
        # retrive information on breakpoints and define new segments (see the function for better explanation)
        # for the desidered type of mutations
        multi_cna <- join_segments(cnaqc_objs = cnaqc_objs, cna_type)
        
        # map the original mutations on the new segments for each segment
        multi_mutations <- lapply(names(multi_cna), function(x) {
            set_elements(cnaqc_obj = cnaqc_objs[[x]], new_cna = multi_cna[[x]], cna_type = cna_type)
        })
        
        names(multi_mutations) = names(multi_cna)
        
        # create a list with new cnaqc objects with the new segmentation and the mutations mapped on them
        multi_input = lapply(names(multi_cna), function(x) {
            init(mutations = multi_mutations[[x]], 
                cna = multi_cna[[x]], 
                purity = cnaqc_objs[[x]]\$purity, 
                sample = x)
        })
        
        names(multi_input) <- lapply(multi_input, function(x) {x\$sample}) %>% unlist()
        
        # assign m_cnaqc class to the output
        class(multi_input) <- "m_cnaqc"
        
        return(multi_input)
    }

    get_segment_info = function(data, chr, sample, new_from, new_to, keep_columns){
        data %>%
            dplyr::filter(sample_id == sample,
                        chr == !!chr,
                        from <= new_from,
                        to >= new_to) %>%
            select(all_of(keep_columns))
    }


    join_segments = function(cnaqc_objs, cna_type){
        # Row binded segments table (with sample specification)
        x = lapply(cnaqc_objs %>% names(), function(x){
            
            CNA(cnaqc_objs[[x]], type = cna_type) %>% 
            dplyr::mutate(sample_id = x)
            
        }) %>%
            do.call(dplyr::bind_rows, .) %>%
            dplyr::select(sample_id, dplyr::everything())
        
        out = lapply(x\$chr %>% unique(), function(chr) {
            new_breakpoints = c(
                x %>%
                    dplyr::filter(chr == !!chr) %>%
                    dplyr::pull(from),
                x %>%
                    dplyr::filter(chr == !!chr) %>%
                    dplyr::pull(to)) %>% 
            unique() %>%
            sort()
            
            #  Separate new breakpoints into segment from and to values
            new_from = new_breakpoints[ !new_breakpoints == dplyr::last(new_breakpoints)] # last element will not be included in the from column 
            new_to = new_breakpoints[!new_breakpoints == dplyr::first(new_breakpoints)] # first element will not be included in the to column 
            
            # iterate on the new breakpoints to subset the cna piled up 
            lapply(new_breakpoints[-1] %>% seq_along(), function(i) {
            
            # iterate for each sample 
            lapply(x\$sample_id %>% unique(), function(s) {
                
                # define which columns must be kept in the new table
                not_wanted = c("from", "to", "length", "size", "segment_id", "chr", "sample_id", "n")
                wanted = setdiff(colnames(x), not_wanted)
                
                
                # get the information for the sample in the new segment 
                tmp = get_segment_info(x, 
                                    chr = chr, 
                                    sample = s, 
                                    new_from = new_from[i], 
                                    new_to = new_to[i], 
                                    keep_columns = wanted)
                
                # do some checking on the result
                if (nrow(tmp) == 0) { # there is no information on copy number on the new segment: insert NA as value of all the columns, except from, to, segment_id and sample_id
                
                tmp_v2 = rep(NA, ncol(tmp))  
                names(tmp_v2) = colnames(tmp)
                
                tmp = tmp_v2 %>% 
                    tibble::as_tibble_row()
                }
                
                # create a tibble with the information on the new breakpoints and include the previously retrieved information 
                tidyr::tibble(chr = chr, 
                            from = new_from[i], 
                            to = new_to[i], 
                            sample_id = s) %>% 
                dplyr::bind_cols(tmp) %>% 
                dplyr::mutate(segment_id = paste(chr, from, to, sep = ":")) 
                
            }) %>% do.call(dplyr::bind_rows, .)
            }) %>% do.call(dplyr::bind_rows, .)
        })  %>% do.call(dplyr::bind_rows, .) # create a unique big tibble with the new segmentation
        
        # remove all the segments that are not correctly shared across samples
        remove_segments = out %>% 
            dplyr::filter(is.na(Major)) %>% 
            dplyr::filter(is.na(minor)) %>% 
            dplyr::pull(segment_id) %>% 
            unique()
        
        all_segments = out %>% 
            pull(segment_id) %>% 
            unique()
        
        keep_segments = setdiff(all_segments, remove_segments)
        
        out = out %>% 
            dplyr::filter(segment_id %in% keep_segments)
        
        # split the cna table by sample id
        out_by_sample = lapply(out\$sample_id %>% unique(), function(x) {
            out %>% 
            dplyr::filter(sample_id == x)
        })
        
        names(out_by_sample) = lapply(out_by_sample, function(x) {x\$sample_id %>% unique}) %>% unlist()
        return(out_by_sample)     
    }

    
    set_elements <- function(cnaqc_obj, new_cna, cna_type) {    
        initial_mutations = Mutations(cnaqc_obj, cna = cna_type)
        new_element = CNAqc:::prepare_input_data(initial_mutations, new_cna, cnaqc_obj\$purity)
        remapped_mut = new_element\$mutations
        return(remapped_mut)
    }

    # Actual code
    result = lapply(strsplit("$rds_list", " ")[[1]], FUN = function(file){
             readRDS(file)
             }) 
    names(result) = samples
    out = multisample_init(result)
    saveRDS(object = out, file = paste0(res_dir, "join_cnaqc.rds"))

    join_mut <- tibble()
    join_cna <- tibble()

    for (s in samples){
        purity <- out[[s]]\$purity

        tmp_mut <- out[[s]]\$mutations
        tmp_mut <- dplyr::bind_cols(tmp_mut, purity = rep(purity, nrow(tmp_mut)))
        join_mut <- dplyr::bind_rows(join_mut, tmp_mut)
        
        tmp_cna <- out[[s]]\$cna
        tmp_cna <- dplyr::bind_cols(tmp_cna, purity = rep(purity, nrow(tmp_cna)))
        join_cna <- dplyr::bind_rows(join_cna, tmp_cna)
    }

    mut_join_table <- join_mut %>%  
        dplyr::mutate(patient_id = rep("$patientID", nrow(.))) %>% 
        dplyr::mutate(normal_cn = rep(2, nrow(.))) %>% 
        tidyr::separate(karyotype, into = c('major_cn', 'minor_cn'), sep = ':') %>% 
        dplyr::mutate(mutation_id = paste(patient_id, chr, POS, alt, sep = ':')) %>% 
        dplyr::select(mutation_id, sample, chr, POS, ref, alt, DP, VAF, NR, NV,normal_cn, major_cn, minor_cn, purity, patient_id, SYMBOL, is_driver, segment_id, -ChromKey, everything()) %>% 
        dplyr::rename(sample_id = sample,
                        pos = POS, 
                        gene = SYMBOL)

    cna_join_table <- join_cna %>% dplyr::mutate(patient_id = rep("$patientID", nrow(.))) 
    
    write.table(mut_join_table,  paste0(res_dir, 'mut_join_table.tsv'), sep = "\t", row.names = F, col.names = T, quote = F)
    write.table(cna_join_table,  paste0(res_dir, 'cna_join_table.tsv'), sep = "\t", row.names = F, col.names = T, quote = F)
    """
}
