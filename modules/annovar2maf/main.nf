process ANNOVAR2MAF {
    publishDir params.publish_dir, mode: 'copy'

    
    input:

      tuple val(datasetID), val(patientID), val(sampleID), path(vcf_File)

    output:

      tuple val(datasetID), val(patientID), val(sampleID), path("$datasetID/$patientID/$sampleID/annovar2maf/*.maf")

    script:

    """
    #!/usr/bin/env Rscript

    
    library(maftools)
    
    dir.create(paste0("$datasetID", "/", "$patientID","/","$sampleID","/annovar2maf"), recursive = TRUE)

    annovar_txt <- read.delim("$vcf_File")
    Tumor_Sample_Barcode = paste0("$patientID","_","$sampleID")
    annovar_txt <- cbind(annovar_txt, Tumor_Sample_Barcode)
   
 
    annovar.hg38_multianno.maf <- maftools::annovarToMaf(annovar = annovar_txt,
                            Center = 'GU',
                            refBuild = "hg38",
                            tsbCol = 'Tumor_Sample_Barcode',
                            table = "refGene",
                            MAFobj = TRUE
                            ) 

    
    save(annovar.hg38_multianno.maf, file=paste0("$datasetID","/","$patientID","/","$sampleID","/annovar2maf/annovar.hg38_multianno.maf"))

    """
}
