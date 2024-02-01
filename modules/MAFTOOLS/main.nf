process MAFTOOLS {
    
    publishDir params.publish_dir, mode: 'copy'

    input:

      tuple val(datasetID), val(patientID), val(sampleID), path(maf_File, stageAs: 'data_vep*.maf')

    output:

      tuple val(datasetID), path("$datasetID/MAFTOOLS/*.pdf"), path("$datasetID/MAFTOOLS/*.rds")

    script:

    """
    #!/usr/bin/env Rscript

    library(maftools)

    dir.create(paste0("$datasetID","/MAFTOOLS"), recursive = TRUE)

    mafs <- lapply(X = strsplit("$maf_File", " ")[[1]], FUN = maftools::read.maf)

    maf_merged = maftools:::merge_mafs(maf = mafs, verbose = TRUE)

    pdf(file = paste0("$datasetID","/MAFTOOLS/maf_summary.pdf"))
    plotmafSummary(maf = maf_merged, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
    dev.off()

    pdf(file = paste0("$datasetID","/MAFTOOLS/oncoplot.pdf"))
    oncoplot(maf = maf_merged, top = 10, removeNonMutated = TRUE)
    dev.off()

    saveRDS(object = maf_merged, file = paste0("$datasetID","/MAFTOOLS/maf_merged.rds"))

    """
}
