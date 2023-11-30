process VCF2MAF {
    publishDir params.publish_dir, mode: 'copy'



    input:

      tuple val(datasetID), val(patientID), val(sampleID), path(vcf_File)

    output:

      tuple val(datasetID), val(patientID), val(sampleID), path("$datasetID/$patientID/$sampleID/VCF2MAF/*.maf")

    script:

    """

    mkdir -p $datasetID/$patientID/$sampleID/VCF2MAF

    gunzip -c $vcf_File > data_vep.vcf


    vcf2maf.pl \\
    --input-vcf data_vep.vcf \\
    --output-maf $datasetID/$patientID/$sampleID/VCF2MAF/data_vep.maf \\
    --tumor-id ${patientID}_${sampleID} \\
    --ref-fasta $params.ref_genome \\
    --vep-data $params.vep_cache_version \\
    --inhibit-vep \\
    --ncbi-build $params.assembly

    """
}
