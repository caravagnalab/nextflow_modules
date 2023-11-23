process MAFTOOLS {
    publishDir params.publish_dir, mode: 'copy'


    input:

      tuple val(patientID), val(sampleID), path(maf)

    output:

      path("maftools/mafplot_vep.pdf"), path("maftools/mafplot_vep.maf")

    script:

    """
    #!/usr/bin/env Rscript

    # Import libraries
    library("maftools")

    #Create output directory
    dir.create(paste0("maftools", recursive = TRUE)

    #Creating maftools plots
    
    maf_vep <- read.maf(maf = "$maf")
    maf_merged = maftools:::merge_mafs(maf=c($maf), verbose = TRUE)

    pdf(file = paste0("maftools/","mafplot_vep.pdf"))
    plotmafSummary(maf = maf_merged, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
    dev.off()

    saveRDS(object = maf_merged, file = paste0("maftools/maf_merged.maf.rds")
    
    """
}
