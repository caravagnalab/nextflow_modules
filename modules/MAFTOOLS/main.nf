//
// Summarize and visualize MAF files
//

process MAFTOOLS {
    
    publishDir params.publish_dir, mode: 'copy'

    input:

      tuple val(datasetID), val(patientID), val(sampleID), path(maf_File, stageAs: 'data_vep*.maf')

    output:

      tuple val(datasetID), path("VariantAnnotation/MAFTOOLS/$datasetID/*.pdf"), path("VariantAnnotation/MAFTOOLS/$datasetID/*.rds")

    script:

    """
    #!/usr/bin/env Rscript

    library(maftools)

    #Create results output directory

    dir.create(paste0("VariantAnnotation/MAFTOOLS/", "$datasetID"), recursive = TRUE)
    
    #Reading maf files
    mafs <- lapply(X = strsplit("$maf_File", " ")[[1]], FUN = maftools::read.maf)

    #Create a maf multisample object, merge multiple maf files.
    #MAF object contains main maf file, summarized data and an oncomatrix which is useful to plot oncoplots.
    maf_merged = maftools:::merge_mafs(maf = mafs, verbose = TRUE)

    #Creating a pdf file for summary output
    pdf(file = paste0("$datasetID","/MAFTOOLS/maf_summary.pdf"))
    
    #Plotting MAF summary   
    plotmafSummary(maf = maf_merged, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
    dev.off()
    
    #Plotting oncoplot
    pdf(file = paste0("$datasetID","/MAFTOOLS/oncoplot.pdf"))
    oncoplot(maf = maf_merged, top = 10, removeNonMutated = TRUE)
    dev.off()
    
    #Saving results object
    saveRDS(object = maf_merged, file = paste0("$datasetID","/MAFTOOLS/maf_merged.rds"))

    """
}
