process MAFTOOLS {
    
    publishDir params.publish_dir, mode: 'copy'

    input:

      tuple val(datasetID), val(patientID), val(sampleID), path(maf_File, stageAs: 'data_vep*.maf')

    output:

      tuple val(datasetID), path("VariantAnnotation/MAFTOOLS/$datasetID/*.pdf"), path("VariantAnnotation/MAFTOOLS/$datasetID/*.rds")

    script:

      def args                              = task.ext.args                                 ?: ''
      def rmOutlier                         = args!='' && args.rmOutlier                    ? "$args.rmOutlier" : ""
      def addStat                           = args!='' && args.addStat                      ? "$args.addStat" : ""
      def dashboard                         = args!='' && args.dashboard                    ? "$args.dashboard" : ""
      def titvRaw                           = args!='' && args.titvRaw                      ? "$args.titvRaw" : ""
      def showBarcodes                      = args!='' && args.showBarcodes                 ? "$args.showBarcodes" : ""     
      def top                               = args!='' && args.top                          ? "$args.top" : "" 
      def minMut                            = args!='' && args.minMut                       ? "$args.minMut" : ""
      def genes                             = args!='' && args.genes                        ? "$args.genes" : ""
      def altered                           = args!='' && args.altered                      ? "$args.altered" : ""
      def removeNonMutated                  = args!='' && args.removeNonMutated             ? "$args.removeNonMutated" : ""


   
    """
    #!/usr/bin/env Rscript

    library(maftools)

    dir.create(paste0("VariantAnnotation/MAFTOOLS/", "$datasetID"), recursive = TRUE)


    mafs <- lapply(X = strsplit("$maf_File", " ")[[1]], FUN = maftools::read.maf)

    maf_merged = maftools:::merge_mafs(maf = mafs, verbose = TRUE)

    pdf(file = paste0("VariantAnnotation/MAFTOOLS/","$datasetID","/maf_summary.pdf"))

    plotmafSummary(maf = maf_merged,
                   rmOutlier = as.logical("$rmOutlier"), 
                   addStat = eval(parse(text="$addStat")), 
                   dashboard = as.logical("$dashboard"), 
                   titvRaw = as.logical("$titvRaw"),
                   showBarcodes = as.logical("$showBarcodes"),
                   top = as.integer("$top"))
    dev.off()

    pdf(file = paste0("$datasetID","/MAFTOOLS/oncoplot.pdf"))

    oncoplot(maf = maf_merged,
             minMut = eval(parse(text="$minMut")),
             genes = eval(parse(text="$genes")),
             altered = as.logical("$altered"),
             top = as.integer("$top"),
             removeNonMutated = as.logical("$removeNonMutated"))

    dev.off()

    saveRDS(object = maf_merged, file = paste0("VariantAnnotation/MAFTOOLS/", "$datasetID","/maf_merged.rds"))

    """
}
