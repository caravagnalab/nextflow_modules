process VCF2MAF {
    publishDir params.publish_dir, mode: 'copy'



    input:

      tuple val(patientID), val(sampleID), path(vcf_File)

    output:

      tuple val(patientID), val(sampleID),  path("$patientID/$sampleID/VCF2MAF/data_vep.maf"), path("combined_vep.maf.gz")

    script:

    """

    mkdir -p $patientID/$sampleID/VCF2MAF

    gunzip -c $vcf_File > data_vep.vcf


    vcf2maf.pl \\
    --input-vcf data_vep.vcf \\
    --output-maf $patientID/$sampleID/VCF2MAF/data_vep.maf \\
    --tumor-id ${patientID}_${sampleID} \\
    --ref-fasta $params.ref_genome \\
    --vep-data $params.vep_cache_version \\
    --inhibit-vep \\
    --ncbi-build $params.assembly

    cat $patientID/$sampleID/VCF2MAF/data_vep.maf $patientID/$sampleID/VCF2MAF/data_vep.maf >> combined_vep.maf
    gzip combined_vep.maf
    """
}
