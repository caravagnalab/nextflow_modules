process PLATYPUS_CALL_VARIANTS {

    publishDir params.publish_dir, mode: 'copy'

    input:

    tuple val(datasetID), val(patientID), val(sampleID), path(vcf_list, stageAs: '*.vcf.gz'), path(tbi_list, stageAs: '*.vcf.gz.tbi') 
    tuple val(datasetID), val(patientID), val(sampleID), path(normal_bam_list, stageAs: 'normal.bam'), path(normal_bai_list, stageAs: 'normal.bam.bai')
    tuple val(datasetID), val(patientID), val(sampleID), path(tumor_bam_list, stageAs: '*.bam'), path(tumor_bai_list, stageAs: '*.bam.bai')

    output:

    tuple val(datasetID), val(patientID), val(sampleID), path("$datasetID/$patientID/Platypus/*.vcf"), emit: vcf

    script:

    """
    #!/bin/bash

    out_dir=$datasetID/$patientID/Platypus
    mkdir -p \$out_dir


    """
}