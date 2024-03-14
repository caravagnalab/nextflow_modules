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


process BCFTOOLS_MPILEUP {
    publishDir params.publish_dir, mode: 'copy'

    input:

    tuple val(datasetID), val(patientID), val(sampleID), path(vcf_list, stageAs: '*.vcf.gz'), path(tbi_list, stageAs: '*.vcf.gz.tbi') 
    tuple val(datasetID), val(patientID), val(sampleID), path(normal_bam_list, stageAs: 'normal.bam'), path(normal_bai_list, stageAs: 'normal.bam.bai')
    tuple val(datasetID), val(patientID), val(sampleID), path(tumor_bam_list, stageAs: '*.bam'), path(tumor_bai_list, stageAs: '*.bam.bai')

    output:

    tuple val(datasetID), val(patientID), val(sampleID), path("$datasetID/$patientID/Platypus/*.vcf"), emit: vcf

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mpileup = save_mpileup ? "| tee ${prefix}.mpileup" : ""
    def bgzip_mpileup = save_mpileup ? "bgzip ${prefix}.mpileup" : ""
    def intervals = intervals ? "-T ${intervals}" : ""
    """
    echo "${meta.id}" > sample_name.list

    bcftools \\
        mpileup \\
        --fasta-ref $fasta \\
        $args \\
        $bam \\
        $intervals \\
        $mpileup \\
        | bcftools call --output-type v $args2 \\
        | bcftools reheader --samples sample_name.list \\
        | bcftools view --output-file ${prefix}.vcf.gz --output-type z $args3

    $bgzip_mpileup

    tabix -p vcf -f ${prefix}.vcf.gz

    bcftools stats ${prefix}.vcf.gz > ${prefix}.bcftools_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
