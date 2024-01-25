process ANNOVAR_ANNOTATE {
    publishDir params.publish_dir, mode: 'copy'


    input:

      tuple val(datasetID), val(patientID), val(sampleID), path(vcf_File)


    output:

      tuple val(datasetID), val(patientID), val(sampleID), path("$datasetID/$patientID/$sampleID/ANNOVAR/annovar.hg38_multianno.txt")

    script:

    """

    mkdir -p $datasetID/$patientID/$sampleID/ANNOVAR
    gzip -dc $vcf_File > vcf_File



    perl $params.db/table_annovar.pl \\
    vcf_File \\
    $params.humandb \\
    -buildver $params.buildver \\
    -out $datasetID/$patientID/$sampleID/ANNOVAR/annovar \\
    -protocol refGene \\
    -operation g \\
    -nastring . \\
    -vcfinput \\
    -polish \\
    -remove \\
    -thread 4

    awk '{print (NR>1?"$patientID/$sampleID":"Tumor_Sample_Barcode") "\t" \$0}' $datasetID/$patientID/$sampleID/ANNOVAR/annovar.hg38_multianno.txt > $datasetID/$patientID/$sampleID/ANNOVAR/annovar.txt
    """
}
