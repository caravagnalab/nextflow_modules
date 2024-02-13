//
// Mutational signature extraction with SparseSignatures
//

process SPARSE_SIGNATURES {
    publishDir params.publish_dir, mode: 'copy'

    input:

      tuple val(datasetID), path()

    output:

      tuple val(datasetID), path("$datasetID/SparseSig/sign_bestConfig.rds"),
                            path("$datasetID/SparseSig/discovered_signatures.pdf"),
    script:

    """
    #!/usr/bin/env Rscript
    
    library("SparseSignatures")
    library("BSgenome.Hsapiens.1000genomes.hs37d5")
    library("tidyverse")
    library("ggplot2")


    res_SparseSig = paste0("SparseSig/")
    dir.create(res_SparseSig, recursive = TRUE)

    #Input dataset : vcf / tsv / csv joint-table multisample
    # sample | chrom | start | end | ref | alt

   #Import the constructed data file
   #data(ssm560_reduced)
   #vcf <- read.delim("vcf_multisample.tsv", sep = "\t")

   #Generate the patient vs mutation matrix from mutation data
   bsg = BSgenome.Hsapiens.1000genomes.hs37d5
   data(mutation_categories)
   imported_data = import.trinucleotides.counts(data=ssm560_reduced,reference=bsg)

   #pl1_signatures <- patients.plot(trinucleotides_counts=imported_data,samples="PD10010a")
   #ggplot2::ggsave(plot = pl1_signatures, filename = paste0(res_SparseSig, "data_signatures.pdf"), width = 12, height = 18, units = 'in', dpi = 200)

   #load background signature from COSMIC
   data(background)

   #Define a set of parameters on which to perform the estimation
   data(patients)
   
   #estimate the initial values of beta
   starting_betas = startingBetaEstimation(x=patients,K=3:10, background_signature=background)

   #explore the search space of values for the LASSO penalty in order to make a good choice
   #test different values to sparsify beta
   lambda_range = lambdaRangeBetaEvaluation(x=patients,K=10,beta=starting_betas[[8,1]],
                                         lambda_values=c(0.05,0.10))


   #Find the optimal number of signatures and sparsity level: rely on cross-validation

   cv_out = nmfLassoCV(x=patients, 
                 K=3:10,
                 background_signature=background,
                 nmf_runs=1,
                 lambda_values_alpha=c(0.01, 0.05, 0.1),
                 lambda_values_beta=c(0.01, 0.05, 0.1),
                 cross_validation_repetitions=20,
                 cross_validation_iterations = 10,
                 num_processes=8,
                 iterations = 10)

  saveRDS(object = cv_out, file = "cv_out.rds")

  #Use the resulting matrix to find the combination of parameters that yields the lowest MSE
  min_ii = which(cv_out == min(cv_out), arr.ind = TRUE)
  min_Lambda = rownames(cv_out)[min_ii[1]]
  min_K = colnames(cv_out)[min_ii[2]]
  cat("Minimum MSE at:", min_Lambda, "and", min_K, "\n")

  #Discovering the signatures within the dataset: NMF Lasso
  #compute the signatures for the best configuration, i.e., K = 4.
  beta = starting_betas_example[["5_signatures","Value"]]
  nmf_Lasso_out = SparseSignatures::nmfLasso(x = patients, K = min_K, background_signature = background, beta = beta)

  saveRDS(object = nmf_Lasso_out, file = "sign_bestConfig.rds")

  #signature visualization
  signatures = nmf_Lasso_out$beta
  pl2_signatures <- signatures.plot(beta=signatures, xlabels=FALSE)
  
  #saving results 
  ggplot2::ggsave(plot = plot_signatures, filename = paste0(res_SparseSig, "discovered_signatures.pdf"), width = 12, height = 18, units = 'in', dpi = 200)

    """
}
