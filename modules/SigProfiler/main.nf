//
// Mutational signature extraction with SigProfilerExtractor
//

process SIG_PROFILER {
    publishDir params.publish_dir, mode: 'copy'

    input:

      tuple val(datasetID), path()

    output:

      tuple val(datasetID), path("$datasetID/SIGPROFILER/results/SBS96/SBS96_selection_plot.pdf"),
      path("$datasetID/SIGPROFILER/results/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/SBS96_Decomposition_Plots.pdf"), 
      path("$datasetID/SIGPROFILER/results/SBS96/Samples.txt")
    
    script:

    """
    #!/usr/bin/env python3
     
    import os
    import shutil
    from SigProfilerExtractor import sigpro as sig

    #create directory
    path = '/orfeo/scratch/cdslab/kdavydzenka/SIGPROFILER/'
    os.mkdir(path)

    #chose the data type that you would like to import: "vcf" or "matrix"

    data = sig.importdata("vcf") 


    #extract the signatures (use reference_genome parameter only in the input is "vcf")

    sig.sigProfilerExtractor("vcf", "results", "/orfeo/scratch/cdslab/kdavydzenka/BRCA/BRCA_vcf", reference_genome="GRCh37", minimum_signatures=1, maximum_signatures=10, nmf_replicates=100)

    #save the output results

    source_dir = "/orfeo/scratch/cdslab/kdavydzenka/results"
    dest_dir = "/orfeo/scratch/cdslab/kdavydzenka/SIGPROFILER/results"
    shutil.copytree(source_dir, dest_dir)

    """
}
