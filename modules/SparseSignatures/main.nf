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
                            path("$datasetID/SparseSig/cv_out.rds")
    script:

      def args                              = task.ext.args                                 ?: ''
      def K                                 = args!='' && args.K                            ? "$args.K" : "2:10"
      def background_signature              = args!='' && args.background_signature         ? "$args.background_signature" : "background"
      def beta                              = args!='' && args.beta                         ? "$args.beta" : "NULL"
      def normalize_counts                  = args!='' && args.normalize_counts             ? "$args.normalize_counts" : "TRUE"
      def nmf_runs                          = args!='' && args.nmf_runs                     ? "$args.nmf_runs" : "10"
      def iterations                        = args!='' && args.iterations                   ? "$args.iterations" : "30"
      def max_iterations_lasso              = args!='' && args.max_iterations_lasso         ? "$args.max_iterations_lasso" : "10000"
      def num_processes                     = args!='' && args.num_processes                ? "$args.num_processes" : "Inf"
      def starting_beta                     = args!='' && args.starting_beta                ? "$args.starting_beta" : "starting_betas"
      def cross_validation_entries          = args!='' && args.cross_validation_entries     ? "$args.cross_validation_entries" : "0.01"
      def cross_validation_iterations       = args!='' && args.cross_validation_iterations  ? "$args.cross_validation_iterations" : "5"
      def cross_validation_repetitions      = args!='' && args.cross_validation_repetitions ? "$args.cross_validation_repetitions" : "50"
      def lambda_alpha_values               = args!='' && args.lambda_alpha_values          ? "$args.lambda_alpha_values" : "0"
      def lambda_beta_values                = args!='' && args.lambda_beta_values           ? "$args.lambda_beta_values" : "c(0.01, 0.05, 0.10, 0.20)"
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
  
    input_data <- mutations_multisample %>%
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
    mut_counts = SparseSignatures::import.trinucleotides.counts(data=input_data, reference=bsg)

    #load a reference SBS5 background signature from COSMIC
    data(background)

    #OR, for the human germline-derived signature
    #data(background2)

    #estimate the initial values of beta
    starting_betas = SparseSignatures::startingBetaEstimation(x = mut_counts,
                                        K= as.integer("$K"), #user defined
                                        background_signature  = "$background_signature")

    
    #Find the optimal number of signatures and sparsity level: rely on cross-validation
    # 1 h per repetition
    #higher number of CV repetitions corresponds to more accurate parameter estimates
    cv_out = SparseSignatures::nmfLassoCV(
                 x = mut_counts,
                 K = as.integer("$K"),  #user defined
                 starting_beta = "$starting_beta",
                 background_signature = "$background_signature", #provided by user; ignored if beta is given
                 normalize_counts = "$normalize_counts", #for algorithm stability
                 nmf_runs = "$nmf_runs", #number of iterations to estimate the length(K) matrices beta in case beta is NULL; ignored if beta is given
                 lambda_values_alpha = "$lambda_values_alpha", #regularization for the exposures α
                 lambda_values_beta = "$lambda_values_beta",
                 cross_validation_entries = "$cross_validation_entries", #cross-validation test size, i.e., the percentage of entries set to zero during NMF and used for validation
                 cross_validation_iterations = "$cross_validation_iterations", #number of randomized restarts of a single cross-validation repetition, in case of poor fits
                 cross_validation_repetitions = "$cross_validation_repetitions", #number of repetitions of the cross-validation procedure
                 iterations = "$iterations", #number of iterations of every single run of NMF LASSO
                 max_iterations_lasso = "$max_iterations_lasso", #number of sub-iterations involved in the sparsification phase, within a full NMF LASSO iteration
                 num_processes = "$num_processes", #number of requested NMF worker subprocesses to spawn. If Inf (an adaptive maximum number is automatically chosen); if NA or NULL, the function is run as a single process
                 seed = "$seed", verbose = "$verbose", log_file = "$log_file)

    saveRDS(object = cv_out, file = paste0(res_SparseSig/, "cv_out.rds")

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
    min_K <- colnames(cv_means_mse)[min_ii[2]] 
    min_K <- substring(min_K,1,1) %>% as.numeric(min_K)


    #Discovering the signatures within the dataset: NMF Lasso
    #compute the signatures for the best configuration.

    nmf_Lasso_out = SparseSignatures::nmfLasso(
                                           x = mut_counts,
                                           K = min_K,
                                           beta = "$beta", #initial value of the signature matrix β. If NULL, it is estimated with a few runs of NMF.
                                           background_signature = "$background_signature", #provided by the user, is ignored if beta is given instead.If NULL, it is estimated through NMF
                                           normalize_counts = "$normalize_counts",
                                           lambda_rate_alpha = $lambda_alpha_values, #sparsity parameter for the exposure values alpha
                                           lambda_rate_beta = min_Lambda_beta,
                                           iterations = "$iterations", #number of iterations in a single NMF LASSO algorithm run
                                           max_iterations_lasso = "$max_iterations_lasso", #number of sub-iterations involved in the sparsification phase, within a full NMF LASSO iteration
                                           seed = "$seed", verbose = "$verbose")

    saveRDS(object = nmf_Lasso_out, file =  paste0(res_SparseSig, "signatures_bestConfig.rds")

    #signature visualization
    signatures = nmf_Lasso_out$beta
    plot_signatures <- SparseSignatures::signatures.plot(beta=signatures, xlabels=FALSE)

    ggplot2::ggsave(plot = plot_signatures, filename = paste0(res_SparseSig, "disc_signatures.pdf"), width = 12, height = 18, units = 'in', dpi = 200)
  
    """  
}
