process JOIN_CNAQC {
  publishDir params.publish_dir, mode: 'copy'

  input:
    
    tuple val(datasetID), val(patientID), val(sampleID), path(rds_list, stageAs: '*.rds')
  
  output:

    tuple val(datasetID), val(patientID), val(sampleID), path("$datasetID/$patientID/join_CNAqc/*.rds"), emit: rds
    //tuple val(datasetID), val(patientID), val(sampleID), path("$datasetID/$patientID/join_CNAqc/mut_join_table.tsv"), emit: mut_tsv
    //tuple val(datasetID), val(patientID), val(sampleID), path("$datasetID/$patientID/join_CNAqc/cna_join_table.tsv"), emit: cna_tsv

  script:
    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    library(CNAqc)
    
    res_dir = paste0("$datasetID", "/", "$patientID", "/join_CNAqc/")
    dir.create(res_dir, recursive = TRUE)
    
    samples = substr("$sampleID", 2, nchar("$sampleID")-1)
    samples = strsplit(samples, ", ")[[1]]

    
    checking_input <- function(cnaqc_objs) {
        if (class(cnaqc_objs) != "list") {
            wrong_class_all = class(cnaqc_objs)
            cli::cli_abort(c("cnaqc_objs must be a list of CNAqc objects", "x" = "{.var cnaqc_objs} is not a list, you supplied a {.cls {class(cnaqc_objs)}} instead"))
        }
        
        if (unique(lapply(cnaqc_objs, class)) != "cnaqc") {
            wrong_class_elements = unique(unique(lapply(cnaqc_objs, class)))
            cli::cli_abort(c("{.var cnaqc_objs} must be a list of CNAqc objects", "x" = "Elements of {.var cnaqc_objs} are of class {.val {wrong_class_elements}} instead"))
        }
        
        if (names(cnaqc_objs) %>% is.null()) {
            cli::cli_abort(c("{.var cnaqc_objs} is not a named list", "x" = "Names of {.var cnaqc_objs} must correspond to the sample_id"))
        }
    }

    multisample_init <- function(cnaqc_objs, 
                             cna_type = "clonal", 
                             QC_filter = TRUE) {
  
        cli::cli_h1("multi_CNAqc - Defining common segments")

        checking_input(cnaqc_objs)
        
        len = length(cnaqc_objs)
        cli::cli_alert_info("Selected CNA type: {.field {cna_type}}")
        cli::cli_alert_info("Found {.field {len}} CNAqc objects:")
        cli::cli_ul(names(cnaqc_objs)) 
        
        cli::cli_h2("Building a {.cls multi_CNAqc} object")
        cat("\n")
        cli::cli_rule("Selecting new segments")
        
        multi_cna <- join_segments(cnaqc_objs = cnaqc_objs, cna_type)
                
        cli::cli_rule("Mapping mutations on new segments")
        
        multi_mutations <- lapply(names(multi_cna), function(x) {
            set_elements(cnaqc_obj = cnaqc_objs[[x]], new_cna_list = multi_cna[[x]], cna_type = cna_type, QC_filter = QC_filter)
        })
        names(multi_mutations) = names(multi_cna)
        
        # multi_input = prepare_input_data_multiple(cnaqc_objs, cna_type)
        
        # create a list with new cnaqc objects with the new segmentation and the mutations mapped on them
        
        cli::cli_h1("Defining the {.cls multi_CNAqc} object")
        
        # multi_input <- list()
        
        multi_input = lapply(names(multi_cna), function(x) {
            
            shared = init(
                mutations = multi_mutations[[x]]\$shared,
                cna = multi_cna[[x]]\$shared,
                purity = cnaqc_objs[[x]]\$purity,
                sample = x)
            private = multi_mutations[[x]]\$private
            original = setdiff(names(cnaqc_objs[[x]]), names(shared))
            other_info = lapply(original, function(o) {
                o = cnaqc_objs[[x]][[o]]
                })
            names(other_info) = original
            
            list(shared_mutations = shared, 
                private_mutations = private, 
                original_additional_info = other_info)
        })
                
        names(multi_input) <- lapply(multi_input, function(x) {x\$shared\$sample}) %>% unlist()
        class(multi_input) <- "m_cnaqc"
        
        if(exists("multi_input", inherits = F) & length(multi_input) == length(cnaqc_objs)) {
            cli::cli_h1("Ended")
            cli::cli_alert_success("{.cls multi_CNAqc} object created including all samples")
        }
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

    join_segments = function(cnaqc_objs, cna_type = 'clonal', QC_filter = TRUE){
        x = lapply(cnaqc_objs %>% names(), function(x){
            CNA(cnaqc_objs[[x]], type = cna_type) %>% 
            dplyr::mutate(sample_id = x)
            
        }) %>%
            do.call(bind_rows, .) %>%
            dplyr::select(sample_id, dplyr::everything())
        
        if(QC_filter == TRUE) {
            x = x %>% 
            filter(QC_PASS == TRUE)
        } 
        
        out = lapply(x\$chr %>% unique(), function(chr) {
            # cli::cli_alert_info("Iterating on {.field {chr}}")
            # cli::cli_h2("{symbol\$info} Iterating on {.field {chr}}")
            
            old_segments = sapply(x\$sample_id %>% unique(), function(s) {
            x %>%
                dplyr::filter(sample_id == s) %>%
                dplyr::filter(chr == !!chr) %>%
                dplyr::pull(segment_id) %>%
                unique() %>%
                length()
            })
            
            cli::cli_alert_info("Number of original segments in individual {.cls CNAqc} objects in {.field {chr}}:")
            cli::cli_ul(paste(names(old_segments), old_segments, sep = " = "))
            cat("\n")
                
            new_breakpoints = c(
            x %>%
                dplyr::filter(chr == !!chr) %>%
                dplyr::pull(from),
            x %>%
                dplyr::filter(chr == !!chr) %>%
                dplyr::pull(to)) %>% 
            unique() %>%
            sort()
            
            cli::cli_alert_info("Found {.field {length(new_breakpoints)}} breakpoints")


            new_from = new_breakpoints[ !new_breakpoints == dplyr::last(new_breakpoints)] 
            new_to = new_breakpoints[!new_breakpoints == dplyr::first(new_breakpoints)]
            
            lapply(new_breakpoints[-1] %>% seq_along(), function(i) {
            
            # iterate for each sample 
            lapply(x\$sample_id %>% unique(), function(s) {
                not_wanted = c("from", "to", "length", "size", "segment_id", "chr", "sample_id", "n")
                wanted = setdiff(colnames(x), not_wanted)
                
                tmp = get_segment_info(x, 
                                    chr = chr, 
                                    sample = s, 
                                    new_from = new_from[i], 
                                    new_to = new_to[i], 
                                    keep_columns = wanted)
                
                if (nrow(tmp) == 0) { 
                    tmp_v2 = rep(NA, ncol(tmp))  
                    names(tmp_v2) = colnames(tmp)
                    
                    tmp = tmp_v2 %>% 
                        tibble::as_tibble_row()
                }
                
                tidyr::tibble(chr = chr, 
                    from = new_from[i], 
                    to = new_to[i], 
                    sample_id = s) %>% 
                dplyr::bind_cols(tmp) %>% 
                dplyr::mutate(segment_id = paste(chr, from, to, sep = ":")) 
                
                }) %>% do.call(bind_rows, .)
            }) %>% do.call(bind_rows, .)
        }) %>% do.call(bind_rows, .)

        remove_segments = out %>% 
            dplyr::filter(is.na(Major)) %>% 
            dplyr::filter(is.na(minor)) %>% 
            dplyr::pull(segment_id) %>% 
            unique()

        cat("\n")
        cli::cli_alert_warning("Found {.val {length(remove_segments)}} not shared segments")
        
        all_segments = out %>%
            pull(segment_id) %>%
            unique()
        
        keep_segments = setdiff(all_segments, remove_segments)
        
        out_shared = out %>% dplyr::filter(segment_id %in% keep_segments)
        
        out_private = out %>%
            dplyr::filter(segment_id %in% remove_segments) %>%
            dplyr::filter(!is.na(minor) & !is.na(Major))
            
        cli::cli_alert_info("Shared segments across samples: {.val {out_shared %>% pull(segment_id) %>% unique() %>%  length()}}")
        
        out_shared_by_sample = lapply(out_shared\$sample_id %>% unique(), function(x) {
            out_shared %>%
            dplyr::filter(sample_id == x)})
        
        out_private_by_sample = lapply(out_private\$sample_id %>% unique(), function(x) {
            if (nrow(out_private) == 0) {
            out_private %>%
                dplyr::add_row(chr = NA)
            } else {
            out_private %>%
                dplyr::filter(sample_id == x)
            }
        })
        
        if(length(out_private_by_sample) == 0) {
            lapply(x\$sample_id %>% unique(), function(l) {
            out_private_by_sample[[l]] = rep(NA, ncol(out_private)) 
            })
        } else {
            names(out_private_by_sample) = lapply(out_private_by_sample, function(x) {x\$sample_id %>% unique}) %>% unlist()
        }
        
        names(out_shared_by_sample) = lapply(out_shared_by_sample, function(x) {x\$sample_id %>% unique}) %>% unlist()
        names(out_shared_by_sample) = lapply(out_shared_by_sample, function(x) {x\$sample_id %>% unique}) %>% unlist()
  
        out_all_segments = lapply(x\$sample_id %>% unique, function(x) {
            list(shared = out_shared_by_sample[[x]], private = out_private_by_sample[[x]])
        })
        
        names(out_all_segments) = lapply(out_all_segments, function(x) {x\$shared\$sample_id %>% unique()}) %>% unlist()
        
        cli::cli_alert_success("Obtained updated CNA table with shared and private segments per sample")  
        return(out_all_segments)
    }

    set_elements <- function(cnaqc_obj, new_cna_list, cna_type, QC_filter) {
        cli::cli_rule(crayon::bgCyan(crayon::white(cnaqc_obj\$sample)))
        
        if(class(new_cna_list) != "list") {
            cli::cli_abort(c("Provided a {.cls {class(new_cna_list)}} object", 
                            "x" = "Input must be a list with the new segments"))}
        
        initial_mutations = Mutations(cnaqc_obj, cna = cna_type)
        
        if(QC_filter == TRUE) {
            initial_mutations = initial_mutations %>% 
            filter(QC_PASS == TRUE)}
        
        cli::cli_alert_info("Found {.val {nrow(initial_mutations)}} mutations in the original {.cls CNAqc} object")
                
        # cli::cli_rule(paste(crayon::cyan("{symbol\$info}"), "Mapping mutations on shared segments"))
        mutations_shared_segments = CNAqc:::prepare_input_data(mutations = initial_mutations, cna = new_cna_list\$shared, tumour_purity = cnaqc_obj\$purity)
        remapped_mut_shared = mutations_shared_segments\$mutations
        
        cat("\n")
        # cli::cli_rule(paste(crayon::cyan("{symbol\$info}"), "Mapping mutations on private segments"))
        
        if(new_cna_list\$private %>% is.na() %>% all() == TRUE) {
            cli::cli_alert_warning("No private segments found")
            remapped_mut = list(shared = remapped_mut_shared, private = NA)
            
        } else {
            mutations_private_segments = CNAqc:::prepare_input_data(mutations = initial_mutations, cna = new_cna_list\$private, tumour_purity = cnaqc_obj\$purity)
            remapped_mut_private = mutations_private_segments\$mutations
            remapped_mut = list(shared = remapped_mut_shared, private = remapped_mut_private)
            cat("\n")
        }
        return(remapped_mut)
    }

    # Actual code
    result = lapply(strsplit("$rds_list", " ")[[1]], FUN = function(file){
             readRDS(file)
             }) 
    names(result) = samples
    out = multisample_init(result)
    saveRDS(object = out, file = paste0(res_dir, "multi_cnaqc.rds"))
    """
}
