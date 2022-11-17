#!/usr/bin/env Rscript

library(sequenza)

chromosomes = 1:24
if (as.logical("FALSE"))
  chromosomes = 1:23

seqzExt <- sequenza.extract(
  file = "mini.seqz.gz",
  chromosome.list = chromosomes,
  normalization.method = "median",
  window = as.numeric("1e5"),
  gamma = as.integer("280"),
  kmin = as.integer("300),
  min.reads.baf = as.integer("50"),
  min.reads = as.integer("50"),
  min.reads.normal = as.integer("15"),
  max.mut.types = as.integer("1")
  )

paraSpace <- sequenza.fit(
  sequenza.extract = seqzExt,
  cellularity = seq(as.integer("1"), as.integer("1"), 0.01),
  ploidy = seq(as.integer("1"), as.integer("1"), 0.1),
  # chromosome.list = chr.fit,
  chromosome.list = chromosomes,
  female = as.logical("FALSE")
  )

sequenza.results(
  sequenza.extract = seqzExt,
  cp.table = paraSpace,
  sample.id = "",
  out.dir = "/mnt/c/Users/Elena/Desktop/tirocinio/repo/nextflow_modules",  # should it be the process folder?
  female = as.logical("FALSE")
  )
