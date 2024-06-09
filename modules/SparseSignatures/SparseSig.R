#!/usr/bin/env Rscript

library("SparseSignatures")
library("tidyverse")
library("ggplot2")

res_SparseSig = paste0("SparseSig_test/")

dir.create(res_SparseSig, recursive = TRUE)

#Input dataset : vcf / tsv / csv joint-table multisample
# sample | chrom | start | end | ref | alt
input_data <- read.delim(file = '/orfeo/LTS/CDSLab/LT_storage/kdavydzenka/nextflow_modules/modules/SparseSignatures/mutationsTable.tsv', sep = '\t', header = TRUE)
#input_data <- input_data[ (input_data[["chr"]] %in% c("chr1")), ]

#Extract input data information
input_data <- input_data[,c("Indiv","chr","from","to","ref","alt")]
input_data <- setNames(input_data, c("sample","chrom","start","end","ref","alt"))
input_data[["end"]] <- input_data[["start"]]
input_data[["chrom"]] <- substring(input_data[["chrom"]],4,5)

#Generate the patient vs mutation count matrix from mutation data
#Install a reference human-genome specification. 
#The user must select, among the available choices, the reference genome consistent with the mutation dataset.

library("BSgenome.Hsapiens.1000genomes.hs37d5")
bsg = BSgenome.Hsapiens.1000genomes.hs37d5::hs37d5
mut_counts = import.trinucleotides.counts(data=input_data, 
					  reference=bsg)

#load a reference SBS5 background signature from COSMIC
data(background)

#OR, for the human germline-derived signature
#data(background2)

#estimate the initial values of beta
starting_betas = startingBetaEstimation(x = mut_counts, 
					K= 2:10, #user defined
					background_signature  = background)

#Find the optimal number of signatures and sparsity level: rely on cross-validation
# 1 h per repetition
#higher number of CV repetitions corresponds to more accurate parameter estimates
cv_out = nmfLassoCV(x = mut_counts,
		 K = 2:10,  #user defined
		 starting_beta = starting_betas,
		 background_signature = background, #provided by user; ignored if beta is given
		 normalize_counts = TRUE, #for algorithm stability
		 nmf_runs = 10, #number of iterations to estimate the length(K) matrices beta in case beta is NULL; ignored if beta is given
		 lambda_values_alpha = 0, #disabling regularization for the exposures α
		 lambda_values_beta = c(0.01, 0.05, 0.1, 0.2),
		 cross_validation_entries = 0.01, #cross-validation test size, i.e., the percentage of entries set to zero during NMF and used for validation
		 cross_validation_iterations = 5, #number of randomized restarts of a single cross-validation repetition, in case of poor fits
		 cross_validation_repetitions = 10, #number of repetitions of the cross-validation procedure
		 iterations = 30, #number of iterations of every single run of NMF LASSO
		 max_iterations_lasso = 10000, #number of sub-iterations involved in the sparsification phase, within a full NMF LASSO iteration
		 num_processes = Inf, #number of requested NMF worker subprocesses to spawn. If Inf (an adaptive maximum number is automatically chosen); if NA or NULL, the function is run as a single process
		 seed = NULL, verbose = TRUE, log_file = "")

saveRDS(object = cv_out, file = "SparseSig_test/cv_out.rds")


#Analyze the mean squared error results averaging over cross-validation repetitions
cv_mses <- cv_out$grid_search_mse[1, , ]
cv_means_mse <- matrix(sapply(cv_mses, FUN = mean),
		       nrow = dim(cv_mses)[1]
)
dimnames(cv_means_mse) <- dimnames(cv_mses)

#Find the combination of parameters that yields the lowest MSE
min_ii <- which(cv_means_mse == min(cv_means_mse), arr.ind = TRUE)
min_Lambda_beta <- rownames(cv_means_mse)[min_ii[1]] 
min_Lambda_beta <- as.numeric(gsub("_Lambda_Beta", "", min_Lambda_beta))
min_K <- colnames(cv_means_mse)[min_ii[2]]
min_K <- as.numeric(gsub("_Signatures", "", min_K))

cat("Minimum MSE at:", min_Lambda_beta, "and", min_K, "\n")

#Discovering the signatures within the dataset: NMF Lasso
#compute the signatures for the best configuration.

nmf_Lasso_out = SparseSignatures::nmfLasso(x = mut_counts, 
					   K = min_K, 
					   beta = NULL, #initial value of the signature matrix β. If NULL, it is estimated with a few runs of NMF.
					   background_signature = background, #provided by the user, is ignored if beta is given instead.If NULL, it is estimated through NMF
					   normalize_counts = TRUE,
					   lambda_rate_alpha = 0, #sparsity parameter for the exposure values alpha
					   lambda_rate_beta = min_Lambda_beta, 
					   iterations = 30, #number of iterations in a single NMF LASSO algorithm run
					   max_iterations_lasso = 10000, #number of sub-iterations involved in the sparsification phase, within a full NMF LASSO iteration
					   seed = NULL, verbose = TRUE)

saveRDS(object = nmf_Lasso_out, file = "SparseSig_test/signatures_bestConfig.rds")

#signature visualization
signatures = nmf_Lasso_out$beta
plot_signatures <- signatures.plot(beta=signatures, xlabels=FALSE)
ggplot2::ggsave(plot = plot_signatures, filename = paste0(res_SparseSig, "disc_signatures.pdf"), width = 12, height = 18, units = 'in', dpi = 200)
                                                                                                                                      
