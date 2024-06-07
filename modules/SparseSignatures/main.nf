//
// Mutational signature extraction with SparseSignatures
//
i
process SPARSE_SIGNATURES {
  publishDir params.publish_dir, mode: 'copy'
  

  input:
    tuple val(datasetID), path(joint_table)

  output:
    tuple val(datasetID), path("signature_deconvolution/SparseSig/$datasetID/cv_out.rds"), emit: sparsesig_cv_rds
                          path("signature_deconvolution/SparseSig/$datasetID/signatures_bestConfig.rds"), emit: sparsesig_bestConfig_rds
                          path("signature_deconvolution/SparseSig/$datasetID/plot_signatures.pdf"), emit: sparsesig_plot_pdf
                            
  script:

    def args = task.ext.args ?: ""
    def K = args!="" && args.K ? "$args.K" : ""
    def background_signature = args!="" && args.background_signature ? "$args.background_signature" : ""
    def beta = args!="" && args.beta ? "$args.beta" : ""
    def normalize_counts = args!="" && args.normalize_counts ? "$args.normalize_counts" : ""
    def nmf_runs = args!="" && args.nmf_runs ? "$args.nmf_runs" : ""
    def iterations = args!="" && args.iterations ? "$args.iterations" : ""
    def max_iterations_lasso = args!="" && args.max_iterations_lasso ? "$args.max_iterations_lasso" : ""
    def num_processes = args!="" && args.num_processes ? "$args.num_processes" : ""
    def cross_validation_entries = args!="" && args.cross_validation_entries ? "$args.cross_validation_entries" : ""
    def cross_validation_repetitions = args!="" && args.cross_validation_repetitions ? "$args.cross_validation_repetitions" : ""
    def cross_validation_iterations = args!="" && args.cross_validation_iterations ? "$args.cross_validation_iterations" : ""
    def lambda_values_alpha = args!="" && args.lambda_values_alpha ? "$args.lambda_values_alpha" : ""
    def lambda_values_beta = args!="" && args.lambda_values_beta ? "$args.lambda_values_beta" : ""
    def lambda_rate_alpha = args!="" && args.lambda_rate_alpha ? "$args.lambda_rate_alpha" : ""
    def verbose = args!="" && args.verbose ? "$args.verbose" : ""
    def seed = args!="" && args.seed ? "$args.seed" : ""

  """
  #!/usr/bin/env Rscript

  library(SparseSignatures)
  library(ggplot2)
  library(stringr)
  library(BSgenome.Hsapiens.1000genomes.hs37d5)
    
  res_dir = paste0("signature_deconvolution/SparseSig/", "$datasetID", "/")
  dir.create(res_dir, recursive = TRUE)

   
  multisample_table = read.delim("$joint_table")
  multisample_table <- multisample_table[ (multisample_table[["chr"]] %in% c("chr1")), ] 
    
  #Extract input data information

  input_data <- multisample_table[,c("Indiv","chr","from","to","ref","alt")]
  input_data <- setNames(input_data, c("sample","chrom","start","end","ref","alt"))
  input_data[["end"]] <- input_data[["start"]]
  input_data[["chrom"]] <- substring(input_data[["chrom"]],4,5)

  #Generate the patient vs mutation count matrix from mutation data
  #Install a reference human-genome specification.
  #The user must select, among the available choices, the reference genome consistent with the mutation dataset.

  
  bsg = BSgenome.Hsapiens.1000genomes.hs37d5::hs37d5
  mut_counts = SparseSignatures::import.trinucleotides.counts(data=input_data, reference=bsg)

  #load a reference SBS5 background signature from COSMIC
  data(background)

  #estimate the initial values of beta
      
  starting_betas = SparseSignatures::startingBetaEstimation(x = mut_counts, 
                                          K = eval(parse(text="$K")),
                                          background_signature  = background)

    
  #Find the optimal number of signatures and sparsity level: rely on cross-validation
  # 1 h per repetition
  #higher number of CV repetitions corresponds to more accurate parameter estimates
      
  cv_out = SparseSignatures::nmfLassoCV(
              x = mut_counts,
              K = eval(parse(text="$K")),  
              starting_beta = starting_betas,
              background_signature = eval(parse(text="$background_signature")), 
              normalize_counts = as.logical("$normalize_counts"), 
              nmf_runs = as.integer("$nmf_runs"), 
              lambda_values_alpha = eval(parse(text="$lambda_values_alpha")), 
              lambda_values_beta = eval(parse(text="$lambda_values_beta")),
              cross_validation_entries = as.numeric("$cross_validation_entries"), 
              cross_validation_iterations = as.integer("$cross_validation_iterations"), 
              cross_validation_repetitions = as.integer("$cross_validation_repetitions"), 
              iterations = as.integer("$iterations"), 
              max_iterations_lasso = as.integer("$max_iterations_lasso"), 
              num_processes = eval(parse(text="$num_processes")), 
              verbose = as.logical("$verbose"),
              seed = as.integer("$seed")
)

  saveRDS(object = cv_out, file = paste0(res_dir, "cv_out.rds"))

  #Analyze the mean squared error results averaging over cross-validation repetitions
 
  cv_mses <- cv_out[["grid_search_mse"]][1,,]
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
  print(min_K)
  print(min_Lambda_beta)

  #Discovering the signatures within the dataset: NMF Lasso
  #compute the signatures for the best configuration.

  nmf_Lasso_out = SparseSignatures::nmfLasso(
                                        x = mut_counts,
                                        K = min_K,
                                        beta = eval(parse(text="$beta")), 
                                        background_signature = eval(parse(text="$background_signature")), 
                                        normalize_counts = as.logical("$normalize_counts"),
                                        lambda_rate_alpha = eval(parse(text="$lambda_rate_alpha")), 
                                        lambda_rate_beta = min_Lambda_beta,
                                        iterations = as.integer("$iterations"), 
                                        max_iterations_lasso = as.integer("$max_iterations_lasso"), 
                                        verbose = as.logical("$verbose")
)

  saveRDS(object = nmf_Lasso_out, file =  paste0(res_dir, "signatures_bestConfig.rds"))

  #signature visualization

  signatures = nmf_Lasso_out[["beta"]]
  plot_signatures <- SparseSignatures::signatures.plot(beta=signatures, xlabels=FALSE)

  ggplot2::ggsave(plot = plot_signatures, filename = paste0(res_dir, "plot_signatures.pdf"), width = 12, height = 18, units = 'in', dpi = 200)

  """ 
}
