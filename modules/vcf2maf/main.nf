//
// Convert VCF annotated file to MAF
//

process VCF2MAF {
    publishDir params.publish_dir, mode: 'copy'

    input:

      tuple val(datasetID), val(patientID), val(sampleID), path(vep_output)

    output:

      tuple val(datasetID), val(patientID), val(sampleID), path("VariantAnnotation/VCF2MAF/$datasetID/$patientID/$sampleID/*.maf"), emit: vcf2maf_out

    script:

    """

    mkdir -p VariantAnnotation/VCF2MAF/$datasetID/$patientID/$sampleID/


    gunzip -c $vep_output > data_vep.vcf


    vcf2maf.pl \\
    --input-vcf data_vep.vcf \\
    --output-maf VariantAnnotation/VCF2MAF/$datasetID/$patientID/$sampleID/data_vep.maf \\
    --tumor-id ${patientID}_${sampleID} \\
    --ref-fasta $params.ref_genome \\
    --vep-data $params.vep_cache_version \\
    --inhibit-vep \\
    --ncbi-build $params.assembly

    """
}
