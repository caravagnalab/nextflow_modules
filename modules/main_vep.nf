process MAFPLOT_VEP {
    publishDir params.publish_dir, mode: 'copy'


    input:

      tuple val(patientID), val(sampleID), path(vcf_File)

    output:

      tuple val(patientID), val(sampleID), path("mafplot_vep.pdf")

    script:

    """
    #!/usr/bin/env Rscript

    # Import libraries
    library("maftools")

    #Create output directory
    dir.create(paste0("$patientID","/","$sampleID"), recursive = TRUE)

    #Creating maftools plots
    mafplot_vep <- read.maf(maf = paste0("/orfeo/LTS/CDSLab/LT_storage/variant_annotation/res_vcf2maf/","$patientID","/","$sampleID","/combined_vep.maf.gz"))


    pdf(file = paste0("$patientID","/","$sampleID","/","mafplot_vep.pdf"))
    plotmafSummary(maf = mafplot_vep, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
    dev.off()
    
    """
}
