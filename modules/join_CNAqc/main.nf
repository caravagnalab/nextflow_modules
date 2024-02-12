process JOIN_CNAQC {
  publishDir params.publish_dir, mode: 'copy'

  input:
    
    tuple val(datasetID), val(patientID), val(sampleID), path(rds_list, stageAs: '*.rds')
  
  output:

    tuple val(datasetID), val(patientID), val(sampleID), path("$datasetID/$patientID/join_CNAqc/*.rds"), emit: rds
    

  script:


    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    library(CNAqc)
    
    res_dir = paste0("$datasetID", "/", "$patientID", "/join_CNAqc/")
    dir.create(res_dir, recursive = TRUE)
    
    samples = substr("$sampleID", 2, nchar("$sampleID")-1)
    samples = strsplit(samples, ", ")[[1]]

    multisample_init <- function(cnaqc_objs) {
    
    #check if you are passing a list of cnaqc objs
    lapply(cnaqc_objs, function(x) {
        stopifnot(inherits(x, "cnaqc"))
    })
    
    output <- list()
    
    output[["joint_segments"]] <- join_segments(cnaqc_objs = cnaqc_objs)
    output[["joint_table"]] <-
        join_mutations(cnaqc_objs = cnaqc_objs, joint_segments = output[["joint_segments"]])
    output[["metadata"]] <- tibble(sample = cnaqc_objs %>% names(),
                                    purity = sapply(cnaqc_objs, function(x)
                                    x\$purity))
    output[["original_data"]] <- cnaqc_objs
    
    return(output)
    
    }

    join_segments = function(cnaqc_objs){
    # Row binded segments table (with sample specification)
    x = lapply(cnaqc_objs %>% names(), function(x){
        cnaqc_objs[[x]]\$cna\$sample = x
        return(cnaqc_objs[[x]]\$cna)
    }) %>%
        do.call(bind_rows, .) %>%
        dplyr::select(sample, dplyr::everything())
    
    out = lapply(x\$chr %>% unique(), function(chr) {
        
        # Chromosome-specific new breakpoints
        new_breakpoints = c(
        x %>%
            dplyr::filter(chr == !!chr) %>%
            dplyr::pull(from),
        x %>%
            dplyr::filter(chr == !!chr) %>%
            dplyr::pull(to)
        ) %>% unique() %>%  sort()
        
        #  Separate new breakpoints into segment from and to values
        new_from = new_breakpoints[-length(new_breakpoints)]
        new_to = new_breakpoints[-1]
        
        lapply(new_from %>% seq_along(), function(i) {
        lapply(x\$sample %>% unique(), function(s) {
            
            # Get sample-specific segment quantities (Major and minor copy numbers, covRatio)
            get_segment_info = function(x, what){
            x %>%
                dplyr::filter(sample == s,
                            chr == !!chr,
                            from <= new_from[i],
                            to >= new_to[i]) %>%
                dplyr::pull(what)
            }
            
            Major = get_segment_info(x, what = "Major")
            
            minor = get_segment_info(x, what = "minor")
            
            #covRatio =get_segment_info(x, what = "covRatio")
            
            #  Return joint segments table
            tibble(
            chr = chr,
            from = new_from[i],
            to = new_to[i],
            sample = s,
            Major = ifelse(length(Major) == 0, NA, Major), #  Replace missing values with NA
            minor = ifelse(length(Major) == 0, NA, minor), #  Replace missing values with NA
            #covRatio = ifelse(length(covRatio) == 0, NA, covRatio) #  Replace missing values with NA
            ) %>%
            mutate(segment_id = paste(chr, from, to, sep = ":"))
        }) %>% do.call(bind_rows, .)
        }) %>% do.call(bind_rows, .)
    })  %>% do.call(bind_rows, .)
    
    return(out)
    }

    join_mutations = function(cnaqc_objs, joint_segments){
    
    # create the join table with all the mutations for each sample 
    x = lapply(cnaqc_objs %>% names(), function(x){
        cnaqc_objs[[x]]\$mutations\$sample = x
        return(cnaqc_objs[[x]]\$mutations)
    }) %>%
        do.call(bind_rows, .) %>%
        dplyr::select(sample, dplyr::everything()) 
    
    # map the mutations on the new joint segments
    
    out = lapply(x\$sample %>% unique(), function(s) {
        which_segments = joint_segments %>% dplyr::filter(sample == s) #get the new segments for each sample
        
        which_mutations = x %>%
        dplyr::filter(sample == s) %>%
        dplyr::mutate(segment_id = NA,
                        karyotype = NA) #initialize segment_id and karyotype columns (or override if they already exist)
        
        # iterate over the new segments and map the mutations on them 
        lapply(1:nrow(which_segments), function(i) {
        mappable = which(
            which_mutations\$chr == which_segments\$chr[i] &
            which_mutations\$from >= which_segments\$from[i] &
            which_mutations\$to <= which_segments\$to[i]
        )
        
        # fill karyotype and segment_id columns with the information for each mapped mutation. If the mutation does not map any
        # segment, it will be flagged as NA
        which_mutations\$karyotype[mappable] = paste0(which_segments\$Major[i], ':', which_segments\$minor[i]) 
        which_mutations\$segment_id[mappable] = which_segments\$segment_id[i]
        
        return(which_mutations)
        }) %>% do.call(bind_rows, .)
    }) %>% do.call(bind_rows, .) %>%
        dplyr::mutate(karyotype = ifelse(karyotype == "NA:NA", NA, karyotype)) 
    
    return(out)
    }

    result = lapply(strsplit("$rds_list", " ")[[1]], FUN = function(file){
             readRDS(file)
             }) 
    names(result) = samples
    out = multisample_init(result)

    saveRDS(object = out, file = paste0(res_dir, "join_cnaqc.rds"))
    """
}
