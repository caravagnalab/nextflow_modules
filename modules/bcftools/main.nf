process BCFTOOLS_SPLIT_VEP {

    publishDir params.publish_dir, mode: 'copy'

    input:

    tuple val(patientID), val(sampleID), path(vcf_File)

    output:

    path("$patientID/$sampleID/VEP/*.tsv")

    // "bcftools +split-vep snv_vep.vcf.gz -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%CSQ\n' -d -A tab -s worst"
    //  
    // colnames = c("chr","from","ref","alt","FILTER","Allele","Consequence","IMPACT","SYMBOL","Gene","Feature_type","Feature","BIOTYPE","EXON","INTRON","HGVSc","HGVSp","cDNA_position","CDS_position","Protein_position","Amino_acids","Codons","Existing_variation","DISTANCE","STRAND","FLAGS","VARIANT_CLASS","SYMBOL_SOURCE","HGNC_ID","CANONICAL","MANE_SELECT","MANE_PLUS_CLINICAL","TSL","APPRIS","CCDS","ENSP","SWISSPROT","TREMBL","UNIPARC","UNIPROT_ISOFORM","GENE_PHENO","SIFT","PolyPhen","DOMAINS","miRNA","HGVS_OFFSET","AF","AFR_AF","AMR_AF","EAS_AF","EUR_AF","SAS_AF","gnomADe_AF","gnomADe_AFR_AF","gnomADe_AMR_AF","gnomADe_ASJ_AF","gnomADe_EAS_AF","gnomADe_FIN_AF","gnomADe_NFE_AF","gnomADe_OTH_AF","gnomADe_SAS_AF","gnomADg_AF","gnomADg_AFR_AF","gnomADg_AMI_AF","gnomADg_AMR_AF","gnomADg_ASJ_AF","gnomADg_EAS_AF","gnomADg_FIN_AF","gnomADg_MID_AF","gnomADg_NFE_AF","gnomADg_OTH_AF","gnomADg_SAS_AF","MAX_AF","MAX_AF_POPS","CLIN_SIG","SOMATIC","PHENO","PUBMED","MOTIF_NAME","MOTIF_POS","HIGH_INF_POS","MOTIF_SCORE_CHANGE","TRANSCRIPTION_FACTORS","HGVSp")

    script:

    """
    
    mkdir -p $patientID/$sampleID/VEP/

    #!/bin/bash
    
    name=\$(echo $vcf_File | cut -d "." -f 1)
    echo \$name
    bcftools +split-vep $vcf_File -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%CSQ\t%IMPACT\t%SYMBOL\t%Gene\t%EXON\t%HGVSc\t%HGVSp\n' -A "tab" -s worst -o $patientID/$sampleID/VEP/\${name}.tsv

    """
}
