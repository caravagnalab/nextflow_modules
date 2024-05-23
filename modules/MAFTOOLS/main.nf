process MAFTOOLS {
    
    publishDir params.publish_dir, mode: 'copy'

    input:

      tuple val(datasetID), val(patientID), val(sampleID), path(maf_File, stageAs: 'data_vep*.maf')

    output:

      tuple val(datasetID), path("$datasetID/MAFTOOLS/*.pdf"), path("$datasetID/MAFTOOLS/*.rds")

    script:

      def args                              = task.ext.args                                 ?: ''
      def rmOutlier                         = args!='' && args.rmOutlier                    ? "$args.K" : "TRUE"
      def addStat                           = args!='' && args.addStat                      ? "$args.addStat" : "median"
      def dashboard                         = args!='' && args.dashboard                    ? "$args.dashboard" : "TRUE"
      def titvRaw                           = args!='' && args.titvRaw                      ? "$args.titvRaw" : "FALSE"
      def showBarcodes                      = args!='' && args.showBarcodes                 ? "$args.showBarcodes" : "FALSE"     
      def top                               = args!='' && args.top                          ? "$args.top" : "10" 
      def minMut                            = args!='' && args.minMiut                      ? "$args.top" : "NULL"
      def genes                             = args!='' && args.genes                        ? "$args.genes" : "NULL"
      def altered                           = args!='' && args.altered                      ? "$args.altered" : "FALSE"
      def removeNonMutated                  = args!='' && args.removeNonMutated             ? "$args.removeNonMutated" : "TRUE"


   
    """
    #!/usr/bin/env Rscript

    library(maftools)

    dir.create(paste0("$datasetID","/MAFTOOLS"), recursive = TRUE)



    mafs <- lapply(X = strsplit("$maf_File", " ")[[1]], FUN = maftools::read.maf)

    maf_merged = maftools:::merge_mafs(maf = mafs, verbose = TRUE)

    pdf(file = paste0("$datasetID","/MAFTOOLS/maf_summary.pdf"))

    plotmafSummary(maf = maf_merged,
                   rmOutlier = "$rmOutlier", 
                   addStat = "$addStat", 
                   dashboard = "$dashboard", 
                   titvRaw = "$titvRaw",
                   showBarcodes = "$showBarcodes"
                   top = "$top")
    dev.off()

    pdf(file = paste0("$datasetID","/MAFTOOLS/oncoplot.pdf"))

    oncoplot(maf = maf_merged,
             minMut = "$minMut",
             genes = "$genes",
             altered = "$altered",
             top = "$top",
             removeNonMutated = "$removeNonMutated")

    dev.off()

    saveRDS(object = maf_merged, file = paste0("$datasetID","/MAFTOOLS/maf_merged.rds"))

    """
}
