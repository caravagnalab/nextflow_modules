#!/usr/bin/env Rscript

library(sequenza)

chromosomes = 1:24
if (False)
  chromosomes = 1:23

seqzExt <- sequenza.extract(
  file = mini.seqz.gz,
  chromosome.list = chromosomes,
  normalization.method = median,
  window = 1E+5,
  gamma = 280,
  kmin = 300,
  min.reads.baf = 50,
  min.reads = 50,
  min.reads.normal = 15,
  max.mut.types = 1
  )

paraSpace <- sequenza.fit(
  sequenza.extract = seqzExt,
  cellularity = seq(1, 1, 0.01),
  ploidy = seq(1, 1, 0.1),
  # chromosome.list = chr.fit,
  chromosome.list = chromosomes,
  female = as.logical(False)
  )

sequenza.results(
  sequenza.extract = seqzExt,
  cp.table = paraSpace,
  sample.id = ,
  out.dir = ,  # should it be the process folder?
  female = as.logical(False)
  )
