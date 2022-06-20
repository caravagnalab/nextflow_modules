// Extract
 
process SEQUENZA_EXTRACT {

    input:
    path seqzFile

    output:
    path XXX

    script:
    """
    #!/usr/bin/env Rscript
    
    library(sequenza)
   
    seqzExt <-
    sequenza.extract(
      file = $seqzFile,
      chromosome.list = chromosomes,
      normalization.method = 'median',
      window = 1e5,
      gamma = 280,
      kmin = 300,
      min.reads.baf = 50,
      min.reads = 50,
      min.reads.normal = 15,
      max.mut.types = 1
  )

   paraSpace <-
    sequenza.fit(
      sequenza.extract = seqzExt,
      cellularity = seq(low_cell, up_cell, 0.01),
      ploidy = seq(low_ploidy, up_ploidy, 0.1),
      chromosome.list = chr.fit,
      female = as.logical(is_female)
  )

   sequenza.results(
     sequenza.extract = seqzExt,
     cp.table = paraSpace,
     sample.id = sample_id,
     out.dir = "results/",
     female = as.logical(is_female)
  )
    """
}

