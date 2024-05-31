process JOIN_POSITIONS {
    publishDir params.publish_dir, mode: 'copy'

    input:
    tuple val(datasetID), val(patientID), val(sampleID), path(rds), path(vcf_pileup)
    tuple val(d2), val(p2), val(s2), path(positions)
 
    output:

    tuple val(datasetID), val(patientID), val(sampleID), path("lifter/mpileup/$datasetID/$patientID/$sampleID/*.rds"), emit: rds

    script:
    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    library(vcfR)

    res_dir = paste0("lifter/mpileup/", "$datasetID", "/", "$patientID", "/", "$sampleID", "/")
    dir.create(res_dir,recursive = T, showWarnings = F)

    all_positions = read.table("$positions", header = T) %>% 
                    dplyr::as_tibble() %>% 
                    dplyr::mutate(id = paste(chr, from, to, sep = ':')) %>% 
                    dplyr::select(id, alt)
    
    
    vcf = readRDS("$rds")
    mutations = vcf[["$sampleID"]]\$mutations

    pileup = vcfR::read.vcfR("$vcf_pileup")
    tb = vcfR::vcfR2tidy(pileup)

    gt_field = tb\$gt %>% 
            tidyr::separate(gt_AD, sep = ',', into = c("NR", "NV")) %>%
            dplyr::rename(DP = gt_DP) %>% 
            dplyr::mutate(
                NR = as.numeric(NR),
                NV = as.numeric(NV),
                VAF = NV/DP) %>% 
            dplyr::mutate(VAF = ifelse(is.na(VAF),0,VAF)) %>% 
            dplyr::rename(sample = Indiv)

    fix_field = tb\$fix %>%
            dplyr::rename(
                chr = CHROM,
                from = POS,
                ref = REF,
                alt = ALT
            ) %>%
            dplyr::rowwise() %>%
            dplyr::mutate(
                from = as.numeric(from) - 1,
                to = from + 1
            ) %>%
            dplyr::ungroup() %>%
            dplyr::select(chr, from, to, ref, alt, AD) 

    if(nrow(fix_field) != nrow(gt_field))
        stop("Mismatch between the VCF fixed fields and the genotypes, will not process this file.")
   
    pileup_mutations = dplyr::bind_cols(fix_field, gt_field)  %>%
                dplyr::select(chr, from, to, ref, NV, DP, VAF, dplyr::everything()) %>%
                dplyr::mutate(id = paste(chr, from, to, sep = ':'))

    bind = inner_join(pileup_mutations, all_positions, by = join_by(id)) %>% select(-id)

    vcf[["$sampleID"]]\$mutations = dplyr::bind_rows(mutations, bind)
    saveRDS(file = paste0(res_dir, "pileup_VCF.rds"), object = vcf)
    """
}
