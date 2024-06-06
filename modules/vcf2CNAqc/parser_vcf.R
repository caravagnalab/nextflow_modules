library(tidyverse)
library(vcfR)


parse_Mutect = function(vcf, sample_id){
    # Transform vcf to tidy 
    tb = vcfR::vcfR2tidy(vcf)

    # Extract gt field and obtain coverage (DP) and variant allele frequency (VAF) fields
    gt_field = tb$gt %>% 
            tidyr::separate(gt_AD, sep = ",", into = c("NR", "NV")) %>%
            dplyr::mutate(
            NR = as.numeric(NR),
            NV = as.numeric(NV),
            DP = NV + NR,
            VAF = NV/DP) %>%  
            dplyr::rename(sample = Indiv)
    
    # Extract sample names
    samples_list = gt_field$sample %>% unique

    # check if VCF is annotated with VEP
    if ("CSQ" %in% tb$meta$ID){
      # VEP specific field extraction
      # Take CSQ field names and split by |
      
      vep_field = tb$meta %>% 
                dplyr::filter(ID == "CSQ") %>% 
                dplyr::select(Description) %>% 
                dplyr::pull() 

    vep_field = strsplit(vep_field, split = "|", fixed = TRUE)[[1]]
    vep_field = vep_field[2:length(vep_field)-1]

    # Tranform the fix field by splittig the CSQ and select the columns needed
    fix_field = tb$fix %>%
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
                        dplyr::select(chr, from, to, ref, alt, CSQ) %>% 
                        tidyr::separate(CSQ, vep_field, sep = "\\|") %>% 
                        dplyr::select(chr, from, to, ref, alt, IMPACT, SYMBOL, Gene) #can add other thing, CSQ, HGSP
    } else {
        # Take from fix field some columns
        fix_field = tb$fix %>%
            dplyr::rename(
            chr = CHROM,
            from = POS,
            ref = REF,
            alt = ALT
            ) %>%
            dplyr::rowwise() %>%
            dplyr::mutate(
            from = as.numeric(from),
            to = from + nchar(alt)
            ) %>%
            dplyr::ungroup() %>%
            dplyr::select(chr, from, to, ref, alt, dplyr::everything(), -ChromKey, -DP) #-DP
    }
    
    # For each sample create the table of mutations 
    calls = lapply(
            samples_list,
            function(s){
            gt_field_s = gt_field %>% dplyr::filter(sample == s)

            if(nrow(fix_field) != nrow(gt_field_s))
                stop("Mismatch between the VCF fixed fields and the genotypes, will not process this file.")

            fits = list()
            fits$sample = s
            fits$mutations = dplyr::bind_cols(fix_field, gt_field_s) %>%
                dplyr::select(chr, from, to, ref, alt, NV, DP, VAF, dplyr::everything())

            fits})

    names(calls) = samples_list
    normal = names(calls)[length(calls)] # last element is always the normal
    calls = calls[c(sample_id, normal)]

    return(calls)
}


parse_Strelka = function(vcf, sample_id){
  calls = NA
  return(calls)
}

parse_Platypus = function(vcf, sample_id){

    tb = vcfR::vcfR2tidy(vcf)

    gt_field = tb$gt %>% dplyr::mutate(
      DP = as.numeric(gt_NR),
      NV = as.numeric(gt_NV),
      VAF = NV/DP
    ) %>% dplyr::rename(sample = Indiv)

    # Extract each sample
    samples_list = gt_field$sample %>% unique

    # check if VCF is annotated with VEP
    if ("CSQ" %in% tb$meta$ID){
      # VEP specific field extraction
      # Take CSQ field names and split by |

      vep_field = tb$meta %>% 
                dplyr::filter(ID == "CSQ") %>% 
                dplyr::select(Description) %>% 
                dplyr::pull() 

    vep_field = strsplit(vep_field, split = "|", fixed = TRUE)[[1]]
    vep_field = vep_field[2:length(vep_field)-1]

    # Tranform the fix field by splittig the CSQ and select the columns needed
    fix_field = tb$fix %>%
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
                        dplyr::select(chr, from, to, ref, alt, CSQ) %>% 
                        tidyr::separate(CSQ, vep_field, sep = "\\|") %>% 
                        dplyr::select(chr, from, to, ref, alt, IMPACT, SYMBOL, Gene) #can add other thing, CSQ, HGSP
    
    } else {
    fix_field = tb$fix %>%
        dplyr::rename(
        chr = CHROM,
        from = POS,
        ref = REF,
        alt = ALT
        ) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
        from = as.numeric(from),
        to = from + nchar(alt)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::select(chr, from, to, ref, alt, dplyr::everything())
    }

    calls = lapply(
        samples_list,
        function(s){
        gt_field_s = gt_field %>% dplyr::filter(sample == s)

        if(nrow(fix_field) != nrow(gt_field_s))
            stop("Mismatch between the VCF fixed fields and the genotypes, will not process this file.")

        fits = list()
        fits$sample = s
        fits$mutations = dplyr::bind_cols(fix_field, gt_field_s) %>%
            dplyr::select(chr, from, to, ref, alt, NV, DP, VAF, dplyr::everything())

        fits
        })

    names(calls) = samples_list

    return(calls)
}