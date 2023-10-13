process JOIN_TABLES {

  publishDir params.publish_dir
  
  input:
    
    tuple val(patientID), val(timepointID), val(sampleID), val(sex), path(seqzFile), path(snv_vcfFile)
  
  output:
  
    path("$patientID/*.rds")

     // colnames = c("chr","from","ref","alt","FILTER","Allele","Consequence","IMPACT","SYMBOL","Gene","Feature_type","Feature","BIOTYPE","EXON","INTRON","HGVSc","HGVSp","cDNA_position","CDS_position","Protein_position","Amino_acids","Codons","Existing_variation","DISTANCE","STRAND","FLAGS","VARIANT_CLASS","SYMBOL_SOURCE","HGNC_ID","CANONICAL","MANE_SELECT","MANE_PLUS_CLINICAL","TSL","APPRIS","CCDS","ENSP","SWISSPROT","TREMBL","UNIPARC","UNIPROT_ISOFORM","GENE_PHENO","SIFT","PolyPhen","DOMAINS","miRNA","HGVS_OFFSET","AF","AFR_AF","AMR_AF","EAS_AF","EUR_AF","SAS_AF","gnomADe_AF","gnomADe_AFR_AF","gnomADe_AMR_AF","gnomADe_ASJ_AF","gnomADe_EAS_AF","gnomADe_FIN_AF","gnomADe_NFE_AF","gnomADe_OTH_AF","gnomADe_SAS_AF","gnomADg_AF","gnomADg_AFR_AF","gnomADg_AMI_AF","gnomADg_AMR_AF","gnomADg_ASJ_AF","gnomADg_EAS_AF","gnomADg_FIN_AF","gnomADg_MID_AF","gnomADg_NFE_AF","gnomADg_OTH_AF","gnomADg_SAS_AF","MAX_AF","MAX_AF_POPS","CLIN_SIG","SOMATIC","PHENO","PUBMED","MOTIF_NAME","MOTIF_POS","HIGH_INF_POS","MOTIF_SCORE_CHANGE","TRANSCRIPTION_FACTORS","HGVSp")
 
  script:
    
    """
    #!/usr/bin/env Rscript
    
    library(tidyverse)
    library(CNAqc)
    library(evoverse)
    
    sample_mutect_SNVs = lapply(snv_vcfFile, function(f){
      evoverse::evoparse_mutect_mutations(f)
    })

    sample_mutect_INDELs = lapply(indel_vcfFile, function(f){
      evoverse::evoparse_mutect_mutations(f)
    })
    
    sample_platypus_mutations = evoverse::evoparse_platypus_mutations(platypus_joint_vcfFile)

    joint_platypus = build_joint_calls(sample_platypus_mutations, sample_sequenza, ref_genome = , one_pass = FALSE)
    joint_mutect = build_joint_calls(sample_mutect_mutations, sample_sequenza, ref_genome = , one_pass = FALSE)
    
    x = genotype_outsample_mutations(joint_mutect, joint_platypus, sample_names)
    
    annotated = add_annotations(x, annotationFiles, files_type = "tsv)

    dir.create("$patientID", recursive = TRUE)

    saveRDS(object = annotated, file = paste0("$patientID","/joint_table.rds"))

    """
}
