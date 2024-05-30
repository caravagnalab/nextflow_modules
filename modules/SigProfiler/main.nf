//
// Mutational signature extraction with SigProfilerExtractor
//

process SIG_PROFILER {
    publishDir params.publish_dir, mode: 'copy'

    input:

      tuple val(datasetID), path($joint_table)

    output:

      tuple val(datasetID), path("$datasetID/SIGPROFILER/results/SBS96/SBS96_selection_plot.pdf"),
      path("$datasetID/SIGPROFILER/results/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/SBS96_Decomposition_Plots.pdf"), 
      path("$datasetID/SIGPROFILER/results/SBS96/Samples.txt")
    
    script:
    
      def args                              = task.ext.args                                 ?: ''
      def reference_genome                  = args!='' && args.reference_genome             ? "$args.reference_genome" : "{GRCh37, GRCh38}"
      def exome                             = args!='' && args.exome                        ? "$args.background_signature" : "False"
      def bed_file                          = args!='' && args.bed_file                     ? "$args.bed_file" : "None"
      def chrom_based                       = args!='' && args.chrom_based                  ? "$args.chrom_based" : "False"
      def plot                              = args!='' && args.plot                         ? "$args.plot" : "False"
      def tsb_stat                          = args!='' && args.tsb_stat                     ? "$args.tsb_stat" : "False"
      def seqInfo                           = args!='' && args.seqInfo                      ? "$args.seqInfo" : "True"
      def cushion                           = args!='' && args.cushion                      ? "$args.cushion" : "100"
      def volume                            = args!='' && args.volume                       ? "$args.volume" : "None"
      def input_type                        = args!='' && args.input_type                   ? "$args.input_type" : "matrix"
      def context_type                      = args!='' && args.context_type                 ? "$args.context_type" : "96,DINUC,ID"
      def minimum_signatures                = args!='' && args.minimum_signatures           ? "$args.minimum_signatures" : "1"
      def maximum_signatures                = args!='' && args.maximum_signatures           ? "$args.maximum_signatures" : "25"
      def nmf_replicates                    = args!='' && args.nmf_replicates               ? "$args.nmf_replicates" : "100"
      def resample                          = args!='' && args.resample                     ? "$args.resample" : "True"
      def seeds                             = args!='' && args.seeds                        ? "$args.seeds" : "random"
      def matrix_normalization              = args!='' && args.matrix_normalization         ? "$args.matrix_normalization" : "gmm"
      def nmf_init                          = args!='' && args.nmf_init                     ? "$args.nmf_init" : "random"
      def min_nmf_iterations                = args!='' && args.min_nmf_iterations           ? "$args.min_nmf_iterations" : "10000"
      def max_nmf_iterations                = args!='' && args.max_nmf_iterations           ? "$args.max_nmf_iterations" : "1000000"
      def nmf_test_conv                     = args!='' && args.nmf_test_conv                ? "$args.nmf_test_conv" : "10000"
      def tolerance                         = args!='' && args.tolerance                    ? "$args.tolerance" : "1e-15"
      def cpu                               = args!='' && args.cpu                          ? "$args.cpu" : "-1"
      def gpu                               = args!='' && args.gpu                          ? "$args.gpu" : "False"
      def batch_size                        = args!='' && args.batch_size                   ? "$args.cpu" : "1"
      def stability                         = args!='' && args.stability                    ? "$args.stability" : "0.8"
      def min_stability                     = args!='' && args.min_stability                ? "$args.min_stability" : "0.2"
      def combined_stability                = args!='' && args.combined_stability           ? "$args.combined_stability" : "1.0"   
      def cosmic_version                    = args!='' && args.cosmic_version               ? "$args.cosmic_version" : "3.1"
      def de_novo_fit_penalty               = args!='' && args.de_novo_fit_penalty          ? "$args.de_novo_fit_penalty" : "0.02"
      def nnls_add_penalty                  = args!='' && args.nnls_add_penalty             ? "$args.nnls_add_penalty" : "0.05"
      def nnls_remove_penalty               = args!='' && args.nnls_remove_penalty          ? "$args.nnls_remove_penalty" : "0.01"
      def initial_remove_penalty            = args!='' && args.initial_remove_penalty       ? "$args.initial_remove_penalty" : "0.05"
      def refit_denovo_signatures           = args!='' && args.refit_denovo_signatures      ? "$args.refit_denovo_signatures" : "True"
      def make_decomposition_plots          = args!='' && args.make_decomposition_plots     ? "$args.make_decomposition_plots" : "True"
      def collapse_to_SBS96                 = args!='' && args.collapse_to_SBS96            ? "$args.collapse_to_SBS96" : "True"
      def get_all_signature_matrices        = args!='' && args.get_all_signature_matrices   ? "$args.get_all_signature_matrices" : "False"
      def export_probabilities              = args!='' && args.export_probabilities         ? "$args.export_probabilities" : "True"
      

    """
    #!/usr/bin/env python3
     
    import os
    import shutil
    from SigProfilerExtractor import sigpro as sig
    from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

    input_path = "$input_path_sigprof/"
    output_path = input_path + "$output_folder_sigprof/"   
    
    if os.path.exists(output_path):
        shutil.rmtree(output_path)
    os.mkdir(output_path)

    #import input data
    input_data = pd.read_csv("$joint_table", sep = '\t')
    
    #input data preprocessing
    def input_processing(data):
        new_columns = {'Project': "$datasetID", 'Genome': '$reference_genome', 'Type': "SOMATIC", 'mut_type': "SNP"}
        df = data.assign(**new_columns)
        df['chr'] = df['chr'].astype(str).str[3:]
        df = df.rename(columns={'Indiv': 'Sample', 'chr': 'chrom', 'from': 'pos_start', 'to': 'pos_end'})
        df["ID"] = df["Sample"]
        df = df.loc[:, ['Project', 'Sample', 'ID', 'Genome', 'mut_type', 'chrom', 'pos_start', 'pos_end', 'ref', 'alt', 'Type']]
        return df
    
    input_data = input_processing(input_data)

    #saving input matrix to txt
    input_data.to_csv('$output_path', sep='\t', index=False, header=True)

    #mutation's counts matrix generation
    input_matrix = matGen.SigProfilerMatrixGeneratorFunc(
            project = "$datasetID", 
            reference_genome = "$reference_genome", 
            path_to_input_file = "$output_path",
            exome = "$exome",
            bed_file = "$bed_file",
            chrom_based = "$chrom_based",
            plot = "$plot",
            tsb_stat = "$tsb_stat",
            seqInfo = "$seqInfo,
            cushion = "$cushion,
            volume = "$volume")

    # Perform model fitting
    sig.sigProfilerExtractor(input_type = "$input_type", 
                             out_put = "results_sigprof", 
                             input_data = "$output_path",  #path to the file
                             context_type = "$contex_type",  
                             exome = "$exome",
                             minimum_signatures = "$minimum_signatures",  
                             maximum_signatures = "$maximum_signatures", 
                             nmf_replicates = "$nmf_replicates", #the number of iteration to be performed to extract each number signature
                             resample = "$resample",
                             matrix_normalization = "$matrix_normalization", 
                             seeds= "$seeds",
                             nmf_init = "$nmf_init", #he initialization algorithm for W and H matrix of NMF
                             min_nmf_iterations = "$min_nmf_iterations", 
                             max_nmf_iterations = "$max_nmf_iterations",
                             nmf_test_conv = "$nmf_test_conv", 
                             nmf_tolerance = "$nmf_tolerance", 
                             cpu = "$cpu",
                             gpu = "$gpu",
                             batch_size = "$batch_size",
                             stability = "$stability",
                             min_stability = "$min_stability",
                             combined_stability = "$combined_stability",
                             cosmic_version = "$cosmic_version",
                             de_novo_fit_penalty = "$de_novo_fit_penalty",
                             nnls_add_penalty = "$nnls_add_penalty",
                             nnls_remove_penalty = "$nnls_remove_penalty",
                             initial_remove_penalty = "$initial_remove_penalty",
                             refit_denovo_signatures = "$refit_denovo_signatures",
                             make_decomposition_plots = "$make_decomposition_plots", 
                             collapse_to_SBS96 = "$collapse_to_SBS96", #SBS288 and SBS1536 Denovo signatures will be mapped to SBS96 reference signatures
                             get_all_signature_matrices = "$get_all_signatures_matrices",
                             export_probabilities = "$export_probabilities"
)
    
    

    #save the output results

    source_dir = ""
    dest_dir = ""
    shutil.copytree(source_dir, dest_dir)

    """
}
