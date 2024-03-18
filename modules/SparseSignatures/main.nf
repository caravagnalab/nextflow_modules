//
// Mutational signature extraction with SparseSignatures
//

process SPARSE_SIGNATURES {
    publishDir params.publish_dir, mode: 'copy'

    input:

      tuple val(datasetID), path(joint_table)

    output:

      tuple val(datasetID), path("$datasetID/SparseSig/signatures_bestConfig.rds"),
                            path("$datasetID/SparseSig/discovered_signatures.pdf"),
                            path("$datasetID/SparseSig/lambda_range_beta.rds")
                            path("$datasetID/SparseSig/lambda_range_alpha.rds")
                            path("$datasetID/SparseSig/cv_out.rds")
    script:

def args                              = task.ext.args                                 ?: ''
def K                                 = args!='' && args.K                            ? "$args.K" : "2:10"
def background_signature              = args!='' && args.background_signature         ? "$args.background_signature" : "background"
def beta                              = args!='' && args.beta                         ? "$args.beta" : "NULL"
def normalize_counts                  = args!='' && args.normalize_counts             ? "$args.normalize_counts" : "TRUE"
def nmf_runs                          = args!='' && args.nmf_runs                     ? "$args.nmf_tuns" : "10"
def iterations                        = args!='' && args.iterations                   ? "$args.iterations" : "30"
def max_iterations_lasso              = args!='' && args.max_iterations_lasso         ? "$args.max_iterations_lasso" : "10000"
def num_processes                     = args!='' && args.num_processes                ? "$args.num_processes" : "Inf"
def starting_beta                     = args!='' && args.starting_beta                ? "$args.starting_beta" : "starting_betas"
def cross_validation_entries          = args!='' && args.cross_validation_entries     ? "$args.cross_validation_entries" : "0.01"
def cross_validation_iterations       = args!='' && args.cross_validation_iterations  ? "$args.cross_validation_iterations" : "5"
def cross_validation_repetitions      = args!='' && args.cross_validation_repetitions ? "$args.cross_validation_repetitions" : "10"
def lambda_rate_alpha                 = args!='' && args.lambda_rate_alpha            ? "$args.lambda_rate_alpha" : "min_Lambda_alpha"
def lambda_rate_beta                  = args!='' && args.lambda_rate_beta             ? "$args.lambda_rate_beta" : "min_Lambda_beta"
def seed                              = args!='' && args.seed                         ? "$args.seed" : "NULL"
def verbose                           = args!='' && args.verbose                      ? "$args.verbose" : "TRUE"
def log_file                          = args!='' && args.log_file                     ? "$args.log_file" : ""
	

    """
#!/usr/bin/env Rscript

library("SparseSignatures")
library("tidyverse")
library("ggplot2")

res_SparseSig = paste0("SparseSig/")
dir.create(res_SparseSig, recursive = TRUE)

#Input dataset : vcf / tsv / csv joint-table multisample
# sample | chrom | start | end | ref | alt
multisample_table <- read.delim(file = '$joint_table', sep = '\t', header = TRUE)

#Extract input data information
input_data <- multisample_table %>%
  mutate(sample = paste(multisample_table$patient_id, "_", multisample_table$sample_id)) %>%
  dplyr::rename(
    sample = sample,
    chrom = chr,
    start = from,
    ref = ref,
    alt = alt) %>%
  mutate(end = start) %>%
  dplyr::select(sample, chrom, start, end, ref, alt) %>%
  as.data.frame()

input_data$chrom=str_sub(input_data$chrom,4,5)

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
starting_betas = SparseSignatures::startingBetaEstimation(x = mut_counts,
                                        K= "$K", #user defined
                                        background_signature  = "$background")

#Determining a valid range for the sparsity parameter
#range of values of the signature sparsity parameter, whose stability for the chosen value of K is to be tested (ensure the convergence of the iterative procedure)
#Decide on a range for K (should be large enough, decided by user, context dependent)
lambda_test_values <- c(0.01, 0.05, 0.1, 0,2)
K_range <- c(2, 3, 4, 5, 6, 7, 8, 9, 10)
lambda_range_beta <- SparseSignatures::lambdaRangeBetaEvaluation(x=mut_counts,
                                         K=K_range, #number of signatures (min=2)
                                         lambda_values = lambda_test_values, # range of values of the signature sparsity parameter
                                         beta = "$beta", #initial value of signature matrix
                                         background_signature = "$background",
                                         normalize_counts = "$normalize_counts", #useful for algorithm stability
                                         nmf_runs = "$nmf_runs", #number of iterations to estimate the length(K) matrices beta (if beta is NULL). Ignored if beta is given
                                         iterations = "$iterations", #number of iterations of every single run of NMF LASSO
                                         max_iterations_lasso = "$max_iterations_lasso", #number of sub-iterations involved in the sparsification phase, within a full NMF LASSO iteration

                                         num_processes = "$num_processes", #number of requested NMF worker subprocesses to spawn. If Inf (an adaptive maximum number is automatically chosen; if NA or NULL (the function is run as a single process)
                                         seed = "$seed", #random number generation; to be set for reproducibility
                                         verbose = "$verbose",
                                         log_file = "$log_file")
#test how the λβ values change in relation to a number of signatures
#values of λβ should guarantees convergence for all K
#Inspect the results manually to verify whether there is a “cutoff” λβ value. If the loglik_progression entries appear to progressively decrease in absolute value, the combination of K and the corresponding lambda-value is feasible.

for (i in 1:length(lambda_test_values)) {
        print(colnames(lambda_range_beta)[[i]])
        print(lambda_range_beta[[i]]$loglik_progression)
}

lambda_range_alpha <- SparseSignatures::lambdaRangeAlphaEvaluation(x=mut_counts,
						 K = K_range,
						 beta = "$beta",
						 background_signature = "$background",
						 normalize_counts = "$normalize_counts",
						 nmf_runs = "nmf_runs",
						 lambda_values = lambda_test_values,
						 iterations = "$iterations",
						 max_iterations_lasso = "$max_iterations_lasso",
						 num_processes = "$num_processes",
						 seed = "$seed",
						 verbose = "$verbose",
						 log_file = "$log_file")
for (i in 1:length(lambda_test_values)) {
        print(colnames(lambda_range_alpha)[[i]])
        print(lambda_range_alpha[[i]]$loglik_progression)
}

#Find the optimal number of signatures and sparsity level: rely on cross-validation
# 1 h per repetition
#higher number of CV repetitions corresponds to more accurate parameter estimates
cv_out = nmfLassoCV(x = mut_counts,
                 K = "$K",  #user defined
                 starting_beta = "$starting_beta",
                 background_signature = "$background", #provided by user; ignored if beta is given
                 normalize_counts = "$normalize_counts", #for algorithm stability
                 nmf_runs = "nmf_runs", #number of iterations to estimate the length(K) matrices beta in case beta is NULL; ignored if beta is given
                 lambda_values_alpha = lambda_test_values, #regularization for the exposures α
                 lambda_values_beta = lambda_test_values,
                 cross_validation_entries = "$cross_validation_entries", #cross-validation test size, i.e., the percentage of entries set to zero during NMF and used for validation
                 cross_validation_iterations = "$cross_validation_iterations", #number of randomized restarts of a single cross-validation repetition, in case of poor fits
                 cross_validation_repetitions = "$cross_validation_repetitions", #number of repetitions of the cross-validation procedure
                 iterations = "$iterations", #number of iterations of every single run of NMF LASSO
                 max_iterations_lasso = "$max_iterations_lasso", #number of sub-iterations involved in the sparsification phase, within a full NMF LASSO iteration
                 num_processes = "$num_processes", #number of requested NMF worker subprocesses to spawn. If Inf (an adaptive maximum number is automatically chosen); if NA or NULL, the function is run as a single process
                 seed = "$seed", verbose = "$verbose", log_file = "$log_file)

saveRDS(object = cv_out, file = paste0(res_SparseSig, "cv_out.rds")
saveRDS(object = lambda_range_beta, file = paste0(res_SparseSig, "lambda_range_beta.rds")
saveRDS(object = lambda_range_alpha, file = paste0(res_SparseSig,"lambda_range_alpha.rds")

#Analyze the mean squared error results averaging over cross-validation repetitions
cv_mses <- cv_out$grid_search_mse[1, , ]
cv_means_mse <- matrix(sapply(cv_mses, FUN = mean),
                       nrow = dim(cv_mses)[1]
)

dimnames(cv_means_mse) <- dimnames(cv_mses)

#Find the combination of parameters that yields the lowest MSE
min_ii <- which(cv_means_mse == min(cv_means_mse), arr.ind = TRUE)
min_Lambda_beta <- rownames(cv_means_mse)[min_ii[1]]
min_Lambda_beta <- substring(min_Lambda_beta,1,3) %>% as.numeric()
min_Lambda_alpha <- min_Lambda_beta
min_K <- colnames(cv_means_mse)[min_ii[2]]
min_K <- substring(min_K,1,1) %>% as.numeric(min_K)
cat("Minimum MSE at:", min_Lambda, "and", min_K, "\n")


#Discovering the signatures within the dataset: NMF Lasso
#compute the signatures for the best configuration.

nmf_Lasso_out = SparseSignatures::nmfLasso(x = mut_counts,
                                           K = min_K,
                                           beta = "$beta", #initial value of the signature matrix β. If NULL, it is estimated with a few runs of NMF.
                                           background_signature = "$background_signature", #provided by the user, is ignored if beta is given instead.If NULL, it is estimated through NMF
                                           normalize_counts = "$normalize_counts",
                                           lambda_rate_alpha = min_Lambda_alpha, #sparsity parameter for the exposure values alpha
                                           lambda_rate_beta = min_Lambda_beta,
                                           iterations = "$iterations", #number of iterations in a single NMF LASSO algorithm run
                                           max_iterations_lasso = "$max_iterations_lasso", #number of sub-iterations involved in the sparsification phase, within a full NMF LASSO iteration
                                           seed = "$seed", verbose = "$verbose")

saveRDS(object = nmf_Lasso_out, file =  paste0(res_SparseSig, "signatures_bestConfig.rds")

#signature visualization
signatures = nmf_Lasso_out$beta
plot_signatures <- signatures.plot(beta=signatures, xlabels=FALSE)

ggplot2::ggsave(plot = plot_signatures, filename = paste0(res_SparseSig, "disc_signatures.pdf"), width = 12, height = 18, units = 'in', dpi = 200)

~                                                                                                                                                        i
