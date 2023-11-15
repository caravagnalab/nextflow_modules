process ANNOVAR_ANNOTATE {
    publishDir params.publish_dir, mode: 'copy'


    input:

      tuple val(patientID), val(sampleID), path(vcf_File)


    output:

      tuple val(patientID), val(sampleID), path("$patientID/$sampleID/ANNOVAR/annovar.hg38_multianno.txt"), path("$patientID/$sampleID/ANNOVAR/annovar.hg38_multianno.vcf.gz")

    script:

    """

    mkdir -p $patientID/$sampleID/ANNOVAR
    gzip -dc $vcf_File > vcf_File



    perl $params.db/table_annovar.pl \\
    vcf_File \\
    $params.humandb \\
    -buildver $params.buildver \\
    -out $patientID/$sampleID/ANNOVAR/annovar \\
    -protocol refGene \\
    -operation g \\
    -nastring . \\
    -vcfinput \\
    -polish \\
    -remove \\
    -thread 4


    bcftools view $patientID/$sampleID/ANNOVAR/annovar.hg38_multianno.vcf -Oz -o  $patientID/$sampleID/ANNOVAR/annovar.hg38_multianno.vcf.gz



    """
}
