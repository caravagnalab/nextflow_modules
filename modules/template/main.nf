process PROCESS_NAME {

    publishDir params.publish_dir, mode: 'copy'

    input:

      tuple val(patientID), val(sampleID), path(input_file)  

    output:

      tuple val(patientID), val(sampleID), path("/path/to/output/file/possibly/including/variables/like/$patientID/$sampleID")

    script:

    """
    #!/usr/bin/env something (e.g. python, Rscript)

    run your_script --input $input_file --output /path/to/output/file/possibly/including/variables/like/$patientID/sample$ID/output.format

    """
}
