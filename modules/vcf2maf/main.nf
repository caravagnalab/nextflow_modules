process VCF2MAF {
    publishDir params.publish_dir, mode: 'copy'

    input:

      tuple val(patientID), val(sampleID), path(vcf_File)

    output:

      tuple val(patientID), val(sampleID), path("$patientID/$sampleID/VEP/*.maf")

    script:

    """

    mkdir -p $patientID/$sampleID/VEP
    
    gunzip -c $vcf_File > data.vcf 

    vcf2maf.pl \\
    --input-vcf data.vcf \\
    --output-maf $patientID/$sampleID/VEP/data_vep.maf \\
    --tumor-id ${patientID}_${sampleID} \\
    --ref-fasta $params.ref_genome \\
    --vep-data $params.vep_cache_version \\
    --ncbi-build $params.assembly

   """
}
