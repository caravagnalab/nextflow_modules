process SNPEFF2MAF {
    publishDir params.publish_dir, mode: 'copy'


    input:

      tuple val(patientID), val(sampleID), path(vcf_File)


    output:

      tuple val(patientID), val(sampleID), path("$patientID/$sampleID/SNPEFF2MAF/snpEff_ann.maf")

    script:

    """

    mkdir -p $patientID/$sampleID/SNPEFF2MAF
    


    java -Xmx8g -jar /orfeo/LTS/CDSLab/LT_storage/variant_annotation/snpEff/snpEff.jar -c /orfeo/LTS/CDSLab/LT_storage/variant_annotation/snpEff/snpEff.config -v GRCh37.75 vcf_File > $patientID/$sampleID/SNPEFF/snpEff_ann.vcf

    """
}
