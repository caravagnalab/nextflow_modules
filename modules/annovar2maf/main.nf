//
// Convert Annovar annotated file to MAF file
//

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
    
    #Create results output directory
    dir.create(paste0("$datasetID", "/", "$patientID","/","$sampleID","/annovar2maf"), recursive = TRUE)
    
    # Converting Annovar output to MAF
 
    annovar_maf <- maftools::annovarToMaf(annovar = "$vcf_File",
                            Center = 'GU',
                            refBuild = "hg38",
                            tsbCol = 'Tumor_Sample_Barcode',
                            table = "refGene",
                            MAFobj = TRUE
                            ) 
    #Saving results
    
    save(annovar_maf, file=paste0("$datasetID","/","$patientID","/","$sampleID","/annovar2maf/annovar.maf"))

    """
}
