process SNPEFF_ANNOTATE {
    publishDir params.publish_dir, mode: 'copy'


    input:

      tuple val(patientID), val(sampleID), path(vcf_File)


    output:

      tuple val(patientID), val(sampleID), path("$patientID/$sampleID/SNPEFF/snpEff_ann.vcf")

    script:

    """

    mkdir -p $patientID/$sampleID/SNPEFF
    gzip -dc $vcf_File > vcf_File


    java -Xmx8g -jar /orfeo/LTS/CDSLab/LT_storage/variant_annotation/snpEff/snpEff.jar -c /orfeo/LTS/CDSLab/LT_storage/variant_annotation/snpEff/snpEff.config -v GRCh37.75 vcf_File > $patientID/$sampleID/SNPEFF/snpEff_ann.vcf

    """
}
