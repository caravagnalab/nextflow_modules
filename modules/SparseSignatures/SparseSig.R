#!/usr/bin/env Rscript

library("SparseSignatures")
library("BSgenome.Hsapiens.1000genomes.hs37d5")
library("tidyverse")
library("ggplot2")


res_SparseSig = paste0("SparseSig/")
dir.create(res_SparseSig, recursive = TRUE)


#Input dataset : vcf / tsv / csv joint-table multisample
# sample | chrom | start | end | ref | alt

#Extract input data information
multisample_table <- read.delim(file = 'mut_join_table.tsv', sep = '\t', header = TRUE)
input_data <- multisample_table %>%
  mutate(sample = paste(multisample_data$patient_id, "_", multisample_data$sample_id)) %>%
  dplyr::rename(
    sample = sample,
    chrom = chr,
    start = from,
    end = to,
    ref = ref,
    alt = alt) %>%
  dplyr::select(sample, chrom, start, end, ref, alt) %>%
  as.data.frame()

#Import the constructed data file
data(ssm560_reduced) #example

#Generate the patient vs mutation count matrix from mutation data
bsg = BSgenome.Hsapiens.1000genomes.hs37d5::hs37d5
data(mutation_categories)
imported_data = import.trinucleotides.counts(data=ssm560_reduced, reference=bsg)

#load a reference background signature from COSMIC
data(background)

#Define a set of parameters on which to perform the estimation
data(patients)

#estimate the initial values of beta
starting_betas = startingBetaEstimation(x = patients, K= 3:10, background_signature  = background)

#Determining a valid range for the sparsity parameter 
#Choice of K is context dependent, choosen by user
#range of values of the signature sparsity parameter, whose stability for the chosen value of K is to be tested
#lambda_test_values <- c(0.01, 0.05, 0.1) 
#lambda_range = lambdaRangeBetaEvaluation(x=patients,
					 #K=5, 
					 #lambda_values = lambda_test_values, 
					 #beta = starting_betas, 
					 #num_processes = Inf) #adaptive maximum number is automatically chosen

#test how the λβ values change in relation to a number of signatures
#values of λβ should guarantees convergence for all K
#Inspect the results manually to verify whether there is a “cutoff” λβ value
#for (i in 1:length(lambda_test_values)) {
	#print(colnames(lambda_range)[[i]])
	#print(lambda_range[[i]]$loglik_progression)

#}

#Find the optimal number of signatures and sparsity level: rely on cross-validation
# 1 h per repetition
cv_out = nmfLassoCV(x = patients,
		 K = 3:10,  #user supplied
		 starting_beta = starting_betas,
		 lambda_values_alpha = 0, #disabling regularization for the exposures α
		 lambda_values_beta = c(0.01, 0.05, 0.1),
		 cross_validation_repetitions = 10,
		 num_processes = Inf) #number of requested NMF worker subprocesses (adaptive maximum number is automatically chosen)

#saveRDS(object = cv_out, file = "cv_out.rds")

#Analyze the mean squared error results averaging over cross-validation repetitions
cv_mses <- cv_out$grid_search_mse[1, , ]
cv_means_mse <- matrix(sapply(cv_mses, FUN = mean),
		       nrow = dim(cv_mses)[1]
)

dimnames(cv_means_mse) <- dimnames(cv_mses)

#Compute the combination with the lowest MSE
min_ii <- which(cv_means_mse == min(cv_means_mse), arr.ind = TRUE)
min_Lambda <- rownames(cv_means_mse)[min_ii[1]] 
min_Lambda <- substring(min_Lambda,1,4) %>% as.numeric()
min_K <- colnames(cv_means_mse)[min_ii[2]]
min_K <- substring(min_K,1,1) %>% as.numeric(min_K)
cat("Minimum MSE at:", min_Lambda, "and", min_K, "\n")


#Discovering the signatures within the dataset: NMF Lasso
#compute the signatures for the best configuration.

nmf_Lasso_out = SparseSignatures::nmfLasso(x = patients, 
					   K = min_K, 
					   beta = NULL, #initial value of the signature matrix β
					   background_signature = background, #provided by the user, is ignored if beta is given instead
					   lambda_rate_alpha = 0,  
					   lambda_rate_beta = min_Lambda, 
					   iterations = 30) #number of iterations in a single NMF LASSO algorithm run

saveRDS(object = nmf_Lasso_out, file = "signatures_bestConfig.rds")

#signature visualization
signatures = nmf_Lasso_out$beta
plot_signatures <- signatures.plot(beta=signatures, xlabels=FALSE)

ggplot2::ggsave(plot = plot_signatures, filename = paste0(res_SparseSig, "disc_signatures.pdf"), width = 12, height = 18, units = 'in', dpi = 200)
                                                                                                                                      
