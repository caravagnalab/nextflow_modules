#!/usr/bin/env nextflow


params.input_path = "/PATH/TO/VCF/DIRECTORY/"
input_path_ch = Channel.of(params.input_path)

params.output_folder = "OUTPUT_FOLDER_NAME"
output_folder_ch = Channel.of(params.output_folder)

process SIGPROFILER {

    input:
    path input_path
    val output_folder
    
    output:

    script:
    """
    #!/usr/bin/env python

    import os
    import glob
    import shutil
    from SigProfilerExtractor import sigpro as sig

    input_path = "$input_path/"
    output_path = input_path + "$output_folder/"

    if os.path.exists(output_path):
        shutil.rmtree(output_path)
    os.mkdir(output_path)

    sig.sigProfilerExtractor(
        input_type = "vcf", 
        output = "$output_folder", 
        input_data = "$input_path", 
        reference_genome="GRCh37", 
        minimum_signatures=1, maximum_signatures=5, 
        nmf_replicates=10, resample = True, batch_size=1, cpu=-1, gpu=False, 
        min_nmf_iterations= 1000, max_nmf_iterations=10000, 
        nmf_test_conv= 1000
    )
    """
}

workflow { 
    SIGPROFILER(input_path_ch, output_folder_ch)
}

