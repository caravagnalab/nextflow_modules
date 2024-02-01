process SIG_PROFILER {
    publishDir params.publish_dir, mode: 'copy'


    input:

      tuple val(datasetID), path(vcf_joint)


    output:

      tuple val(datasetID), path("$datasetID/SIGPROFILER/results/")

    script:

    """
    #!/usr/bin/env python3
     
    import os
    import shutil
    from SigProfilerExtractor import sigpro as sig

    #create directory
    path = '/orfeo/scratch/cdslab/kdavydzenka/SIGPROFILER/'
    os.mkdir(path)

    #chose the data type that you would like to import

    data = sig.importdata("vcf") 


    #extract the signatures
    sig.sigProfilerExtractor("vcf", "results", "/orfeo/scratch/cdslab/kdavydzenka/BRCA/BRCA_vcf", reference_genome="GRCh37", 
    minimum_signatures=1, maximum_signatures=10, nmf_replicates=100)

    #save the output
    source_dir = "/orfeo/scratch/cdslab/kdavydzenka/results"
    dest_dir = "/orfeo/scratch/cdslab/kdavydzenka/SIGPROFILER/results"
    shutil.copytree(source_dir, dest_dir)
    """
}
