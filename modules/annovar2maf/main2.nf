process ANNOVAR2MAF {
    publishDir params.publish_dir, mode: 'copy'


    input:

      tuple val(patientID), val(sampleID), path(vcf_File)
      
    output:

      tuple val(patientID), val(sampleID), path("$patientID/$sampleID/annovar2maf/annovar.hg38_multianno.maf.rds")
      

    script:

    """
    #!/usr/bin/env Rscript


    library(maftools)

    annovar.annotated = paste0("/orfeo/LTS/CDSLab/LT_storage/variant_annotation/results_new/res_annovar/","$patientID/$sampleID/ANNOVAR/annovar.hg38_multianno.txt")

    annovar.hg38_multianno.maf = maftools::annovarToMaf(annovar = annovar.annotated,
                            Center = 'GU',
                            refBuild = "hg38",
                            tsbCol = 'Tumor_Sample_Barcode',
                            table = "refGene",
                            MAFobj = TRUE
                            )

    dir.create(paste0("$patientID","/","$sampleID","/annovar2maf"), recursive = TRUE)

    saveRDS(object = annovar.hg38_multianno.maf, file = paste0("/orfeo/LTS/CDSLab/LT_storage/variant_annotation/","$patientID","/","$sampleID","/annovar2maf/annovar.hg38_multianno.maf.rds"))

    """
}
