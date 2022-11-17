setwd("C:/Users/Elena/Desktop/tirocinio/repo/nextflow_modules")

library(sequenza)
library(tidyverse)

data("example.seqz")

example.seqz %>% class


chromosomes = paste0("", 1:24)
if (as.logical("FALSE"))
  chromosomes = paste0("", 1:23)


write.table(x=example.seqz, file = "./example_data/example.seqz", row.names = F, col.names = T, sep="\t")
pr = read.seqz("./example_data/example.seqz")


seqzExt <- sequenza.extract(file = "./example_data/example.seqz",
                            chromosome.list = chromosomes,
                            # normalization.method = "$norm_method",
                            # window = as.numeric("$window"),
                            # gamma = as.integer("$gamma"),
                            # kmin = as.integer("$kmin"),
                            # min.reads.baf = as.integer("$min_reads_baf"),
                            # min.reads = as.integer("$min_reads"),
                            # min.reads.normal = as.integer("$min_reads_normal"),
                            # max.mut.types = as.integer("$max_mut_types")
                            )

saveRDS(seqzExt, "seqz_prova.Rds")



paraSpace <- sequenza.fit(
  sequenza.extract = seqzExt,
  # cellularity = seq(as.integer("$low_cell"), as.integer("$up_cell"), 0.01),
  # ploidy = seq(as.integer("$low_ploidy"), as.integer("$up_ploidy"), 0.1),
  chromosome.list = chromosomes,
  # female = as.logical("$is_female")
)


sequenza.results(
  sequenza.extract = seqzExt,
  cp.table = paraSpace,
  # sample.id = "$sample_id",
  # out.dir = "$out_dir",  # should it be the process folder?
  # female = as.logical("$is_female")
)



paraSpace_nf = readRDS("./paraSpace_run.Rds")
