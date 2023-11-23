process ANNOVAR2MAF {
    publishDir params.publish_dir, mode: 'copy'

    
    input:

      tuple val(patientID), val(sampleID), path(vcf_File)

    output:

      tuple val(patientID), val(sampleID), path("$patientID/$sampleID/annovar2maf/annovar.hg38_multianno.maf")

    script:

    """
    #!/usr/bin/env Rscript

    
    library(maftools)
    
    dir.create(paste0("$patientID","/","$sampleID","/annovar2maf"), recursive = TRUE)

    annovar_txt <- read.delim(paste0("/orfeo/LTS/CDSLab/LT_storage/variant_annotation/results/","$patientID","/ANNOVAR/annovar.hg38_multianno.txt"), header=TRUE)
    Tumor_Sample_Barcode = paste0("$patientID","_","$sampleID")
    annovar_txt <- cbind(annovar_txt, Tumor_Sample_Barcode)
    write.table(annovar_txt,file = paste0("/orfeo/LTS/CDSLab/LT_storage/variant_annotation/results/","$patientID","/","$sampleID","/ANNOVAR/annovar.hg38_multianno.txt"),sep="\t",row.names=FALSE)

    
    annovar.hg38_multianno.maf <- maftools::annovarToMaf(annovar = paste0("/orfeo/LTS/CDSLab/LT_storage/variant_annotation/results/","$patientID","/","$sampleID","/ANNOVAR/annovar.hg38_multianno.txt"),
                            Center = 'GU',
                            refBuild = "hg38",
                            tsbCol = 'Tumor_Sample_Barcode',
                            table = "refGene",
                            MAFobj = TRUE
                            ) 

    
    save(annovar.hg38_multianno.maf, file=paste0("$patientID","/","sampleID","/annovar2maf/annovar.hg38_multianno.maf"))

    """
}
