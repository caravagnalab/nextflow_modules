process VCF_PROCESSING {
    publishDir params.publish_dir, mode: 'copy'

    input:
     tuple val(datasetID), val(patientID), val(sampleID), path(vcfFile)

    output:
     tuple val(datasetID), val(patientID), val(sampleID), path("$datasetID/$patientID/$sampleID/vcf2CNAqc/*.rds"), emit: rds 

    script:

    """
    #!/usr/bin/env Rscript 
    
    library(tidyverse)
    library(vcfR)
    
    res_dir = paste0("$datasetID", "/", "$patientID", "/", "$sampleID", "/vcf2CNAqc/")
    dir.create(res_dir, recursive = TRUE)

    annotation = "VEP"  
    if (annotation == "VEP"){
        vcf = vcfR::read.vcfR("$vcfFile")
        tb = vcfR::vcfR2tidy(vcf)

        gt_field = tb\$gt %>% 
            tidyr::separate(gt_AD, sep = ',', into = c("NR", "NV")) %>%
            dplyr::mutate(
                NR = as.numeric(NR),
                NV = as.numeric(NV),
                DP = NV + NR,
                VAF = NV/DP) %>% 
            dplyr::rename(sample = Indiv)

        samples_list = gt_field\$sample %>% unique

        vep_field = tb\$meta %>% 
                        dplyr::filter(ID == 'CSQ') %>% 
                        dplyr::select(Description) %>% 
                        pull()
        vep_field = strsplit(vep_field, split="|", fixed=TRUE)[[1]]
        vep_field = vep_field[2:length(vep_field)-1]
        
        fix_field = tb\$fix %>%
                    dplyr::rename(
                        chr = CHROM,
                        from = POS,
                        ref = REF,
                        alt = ALT) %>%
                    dplyr::rowwise() %>%
                    dplyr::mutate(
                        from = as.numeric(from),
                        to = from + nchar(alt)) %>%
                    dplyr::ungroup() %>%
                    dplyr::select(chr, from, to, ref, alt, CSQ, ChromKey, ID, FILTER, DP, NCount) %>% 
                    tidyr::separate(CSQ, vep_field, sep = "\\\\|") %>% 
                    dplyr::select(chr, from, to, ref, alt, FILTER, IMPACT, SYMBOL, Gene) #can add other thing, CSQ, HGSP

        calls = lapply(
            samples_list,
            function(s){
            gt_field_s = gt_field %>% dplyr::filter(sample == s)

            if(nrow(fix_field) != nrow(gt_field_s))
                stop("Mismatch between the VCF fixed fields and the genotypes, will not process this file.")

            fits = list()
            fits\$sample = s
            fits\$caller = "VEP-Mutect"
            fits\$mutations = dplyr::bind_cols(fix_field, gt_field_s) %>%
                            dplyr::select(chr, from, to, ref, alt, NV, DP, VAF, dplyr::everything())
            fits
            })
        names(calls) = samples_list
        
        } else if (annotation == "ANNOVAR"){

            vcf = read.table("$vcfFile", 
                  fill = TRUE, 
                  sep = "\\t", 
                  header = T)

            vcf <- vcf %>% 
                    dplyr::rename(chr = Chr,
                                  from = Start,
                                  to = End,
                                  ref = Ref,
                                  alt = Alt) %>% 
                    dplyr::select(-Otherinfo1, -Otherinfo2, -Otherinfo3, -Otherinfo4, -Otherinfo5, -Otherinfo6, -Otherinfo7, -Otherinfo8, -Otherinfo9) %>% 
                    dplyr::rename(filter = Otherinfo10)
    
            T_res <- apply(vcf, 1,  FUN = function(x){
              values <- strsplit(x['Otherinfo13'][[1]], split = ':') %>% unlist()
              names <- strsplit(x['Otherinfo12'][[1]], split = ':') %>% unlist()
              tidyr::as_tibble(matrix(values, nrow = 1)) %>% 
                `colnames<-`(names) %>% 
                tidyr::separate(AD, sep = ',', into = c('NR', 'NV')) %>% 
                dplyr::rename(VAF = AF)
              }) %>% do.call("bind_rows", .)
            
            N_res <- apply(vcf, 1,  FUN = function(x){
              values <- strsplit(x['Otherinfo14'][[1]], split = ':') %>% unlist()
              names <- strsplit(x['Otherinfo12'][[1]], split = ':') %>% unlist()
              tidyr::as_tibble(matrix(values, nrow = 1)) %>% 
                `colnames<-`(names) %>% 
                tidyr::separate(AD, sep = ',', into = c('NR', 'NV')) %>% 
                dplyr::rename(VAF = AF)
            }) %>% do.call("bind_rows", .)
            
            vcf <- vcf %>% dplyr::select(-Otherinfo12, -Otherinfo13, -Otherinfo14)
            T_res <- bind_cols(vcf, T_res)
            N_res <- bind_cols(vcf, T_res)
            
            tumor <- list()
            tumor$sample <- 'Tumor'
            tumor$caller <- 'ANNOVAR-Mutect'
            tumor$mutations <- T_res
            
            normal <- list()
            normal$sample <- 'Normal'
            normal$caller <- 'ANNOVAR-Mutect'
            normal$mutations <- N_res
            
            calls <- list(tumor, normal)
            }

    saveRDS(object = calls, file = paste0(res_dir, "VCF.rds"))
    """
}
