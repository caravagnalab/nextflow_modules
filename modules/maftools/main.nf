process MAFTOOLS {
    publishDir params.publish_dir, mode: 'copy'


    input:

      tuple val(datasetID), val(patientID), val(sampleID), path(maf_File, stageAs: 'data_vep*.maf')

    output:

      tuple val(datasetID), path("$datasetID/maftools/*.pdf"), path("$datasetID/maftools/*.rds")

    script:

    """
    #!/usr/bin/env Rscript


    library(maftools)
    
    dir.create(paste0("$datasetID","/maftools"), recursive = TRUE)

    mafs <- lapply(X = "$maf_File", FUN = maftools::read.maf)
    
    maf_merged = maftools:::merge_mafs(maf = mafs, verbose = TRUE)

    saveRDS(object = maf_merged, file = paste0("$datasetID","/",maftools/maf_merged.rds")

    pdf(file = paste0("$datasetID","/","maftools/mafplot_vep.pdf"))
    plotmafSummary(maf = maf_merged, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
    dev.off()

    """
}
