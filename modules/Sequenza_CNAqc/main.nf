process SEQUENZA_CNAqc {

  publishDir params.publish_dir
  
  input:
    
    tuple val(patientID), val(timepointID), val(sampleID), val(sex), path(seqzFile), path(snv_vcfFile)
  
  output:
  
    path("$patientID/$timepointID/$sampleID/*.rds")
  
  script:
  
    def args = task.ext.args ?: ''
    def norm_method = args!='' && args.norm_method ? "$args.norm_method" : "median"
    def window = args!='' && args.window ? "$args.window" : "1e5"  // window
    def gamma = args!='' && args.gamma ? "$args.gamma" : "280"  // gamma
    def kmin = args!='' && args.kmin ? "$args.kmin" : "300"  // kmin
    def min_reads_baf = args!='' && args.min_reads_baf ? "$args.min_reads_baf" : "50"
    def min_reads = args!='' && args.min_reads ?  "$args.min_reads" : "50"
    def min_reads_normal = args!='' && args.min_reads_normal ? "$args.min_reads_normal": "15"
    def max_mut_types = args!='' && args.max_mut_types ? "$args.max_mut_types" : "1"
    
    def low_cell = args!='' && args.low_cell ? "$args.low_cell" : "0.95"
    def up_cell = args!='' && args.up_cell ? "$args.up_cell" : "1.0"
    def low_ploidy = args!='' && args.low_ploidy ? "$args.low_ploidy" : "1.8"
    def up_ploidy = args!='' && args.up_ploidy ? "$args.up_ploidy" : "5.4"
    def delta_cellularity = args!='' && args.delta_cellularity ? "$args.delta_cellularity" : "0.05"
    def delta_ploidy = args!='' && args.delta_ploidy ? "$args.delta_ploidy" : "0.25"
    
    """
    #!/usr/bin/env Rscript

    Sys.setenv("VROOM_CONNECTION_SIZE"=99999999)
    
    library(tidyverse)
    library(sequenza)
    library(CNAqc)
    library(evoverse)
    
    SNV = evoverse::evoparse_mutect_mutations("$snv_vcfFile")
   
    SNV = SNV[[1]] 
    SNV\$mutations = SNV\$mutations %>% dplyr::select(chr, from, to, ref, alt, NV, DP, VAF)

    out = CNAqc::Sequenza_CNAqc(
      sample_id = paste0("$patientID"),
      seqz_file = "$seqzFile",
      mutations = SNV\$mutations,
      sex = "$sex",
      cellularity = c(as.integer("$low_cell"), as.integer("$up_cell")),
      ploidy = c(as.integer("$low_ploidy"), as.integer("$up_ploidy")),
      reference = "$params.assembly",
      normalization.method = "$norm_method",
      window = as.numeric("$window"),
      gamma = as.integer("$gamma"),
      kmin = as.integer("$kmin"),
      min.reads.baf = as.integer("$min_reads_baf"),
      min.reads = as.integer("$min_reads"),
      min.reads.normal = as.integer("$min_reads_normal"),
      max.mut.types = as.integer("$max_mut_types"),
      delta_cellularity = as.numeric("$delta_cellularity"),
      delta_ploidy = as.numeric("$delta_ploidy")
    )
    
    dir.create(paste0("$patientID","/","$timepointID","/","$sampleID"), recursive = TRUE)
    saveRDS(object = out, file = paste0("$patientID","/","$timepointID","/","$sampleID","/pipeline.rds"))
    """
}
