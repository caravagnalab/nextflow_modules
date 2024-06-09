process BCFTOOLS_MPILEUP {
    publishDir params.publish_dir
    //mode: 'copy'

    input:

    tuple val(datasetID), val(patientID), val(sampleID), path(bed) 
    tuple val(datasetID), val(patientID), val(sampleID), path(tumor_bam), path(tumor_bai) 

    output:

    tuple val(datasetID), val(patientID), val(sampleID), path("Lifter/mpileup/$datasetID/$patientID/$sampleID/*.vcf"), emit: vcf

    script:
    """
    
    out_dir=lifter/mpileup/$datasetID/$patientID/$sampleID/
    mkdir -p \$out_dir

    bcftools mpileup $tumor_bam \\
        --fasta-ref $params.ref_genome \\
        --regions-file $bed \\
        --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR | bcftools call -Ov -m -o \$out_dir/pileup.vcf

    """
}
