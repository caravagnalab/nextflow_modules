//
// Mutational signature extraction with SparseSignatures
//

process SPARSE_SIGNATURES {
  publishDir params.publish_dir, mode: 'copy'

  input:
    tuple val(datasetID), path(joint_table)

  output:
    tuple val(datasetID), path("$datasetID/SparseSig/*.rds"), emit: rds
                          path("$datasetID/SparseSig/*.pdf"), emit: pdf
                            
  script:

    def args                              = task.ext.args                                 ?: ''
    def K                                 = args!='' && args.K                            ?  "$args.K" : "2:10"
    def background_signature              = args!='' && args.background_signature         ?  "$args.background_signature" : "NULL"
    def beta                              = args!='' && args.beta                         ?  "$args.beta" : "NULL"
    def normalize_counts                  = args!='' && args.normalize_counts             ?  "$args.normalize_counts" : "TRUE"
    def nmf_runs                          = args!='' && args.nmf_runs                     ?  "$args.nmf_runs" : "10"
    def iterations                        = args!='' && args.iterations                   ?  "$args.iterations" : "30"
    def max_iterations_lasso              = args!='' && args.max_iterations_lasso         ?  "$args.max_iterations_lasso" : "10000"
    def num_processes                     = args!='' && args.num_processes                ?  "$args.num_processes" : "Inf"
    def cross_validation_entries          = args!='' && args.cross_validation_entries     ?  "$args.cross_validation_entries" : "0.01"
    def cross_validation_iterations       = args!='' && args.cross_validation_iterations  ?  "$args.cross_validation_iterations" : "5"
    def cross_validation_repetitions      = args!='' && args.cross_validation_repetitions ?  "$args.cross_validation_repetitions" : "10"
    def lambda_values_alpha               = args!='' && args.lambda_values_alpha          ?  "$args.lambda_values_alpha" : "0"
    def lambda_values_beta                = args!='' && args.lambda__values_beta          ?  "$args.lambda_values_beta" : "c(0.01, 0.05, 0.10, 0.20)"
    def verbose                           = args!='' && args.verbose                      ?  "$args.verbose" : "TRUE"

  """
  #!/usr/bin/env Rscript

  library(SparseSignatures)
  library(ggplot2)
  library(stringr)
  library(BSgenome.Hsapiens.1000genomes.hs37d5)
    
  res_SparseSig = paste0("SparseSig/")
  dir.create(res_SparseSig, recursive = TRUE)

   
  multisample_table = read.delim("$joint_table")
   
    
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
              background_signature = "$background_signature", 
              normalize_counts = as.logical("$normalize_counts"), 
              nmf_runs = as.integer("$nmf_runs"), 
              lambda_values_alpha = as.numeric("$lambda_values_alpha"), 
              lambda_values_beta = as.numeric("$lambda_values_beta"),
              cross_validation_entries = as.numeric("$cross_validation_entries"), 
              cross_validation_iterations = as.integer("$cross_validation_iterations"), 
              cross_validation_repetitions = as.integer("$cross_validation_repetitions"), 
              iterations = as.integer("$iterations"), 
              max_iterations_lasso = as.integer("$max_iterations_lasso"), 
              num_processes = "$num_processes", 
              verbose = as.logical($verbose))

  saveRDS(object = cv_out, file = paste0(res_SparseSig, "cv_out.rds"))

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


  #Discovering the signatures within the dataset: NMF Lasso
  #compute the signatures for the best configuration.

  nmf_Lasso_out = SparseSignatures::nmfLasso(
                                        x = mut_counts,
                                        K = min_K,
                                        beta = "$beta", 
                                        background_signature = "$background_signature", 
                                        normalize_counts = as.logical("$normalize_counts"),
                                        lambda_rate_alpha = as.numeric("$lambda_values_alpha"), 
                                        lambda_rate_beta = min_Lambda_beta,
                                        iterations = as.integer("$iterations"), 
                                        max_iterations_lasso = as.integer("$max_iterations_lasso"), 
                                        verbose = as.logical($verbose))

  saveRDS(object = nmf_Lasso_out, file =  paste0(res_SparseSig, "signatures_bestConfig.rds"))

  #signature visualization

  signatures = nmf_Lasso_out[["beta"]]
  plot_signatures <- SparseSignatures::signatures.plot(beta=signatures, xlabels=FALSE)

  ggplot2::ggsave(plot = plot_signatures, filename = paste0(res_SparseSig, "disc_signatures.pdf"), width = 12, height = 18, units = 'in', dpi = 200)

  """ 
}
