//
// Variant Annotation using VEP
//

process VEP_ANNOTATE {
    publishDir params.publish_dir, mode: 'copy'

    input:

      tuple val(datasetID), val(patientID), val(sampleID), path(vcf_File), path(tbi_File)

    output:

      tuple val(datasetID), val(patientID), val(sampleID), path("VariantAnnotation/VEP/$datasetID/$patientID/$sampleID/*.vcf.gz"), emit: vep_output 

    script:

    """

    mkdir -p VariantAnnotation/VEP/$datasetID/$patientID/$sampleID/

    vep --input_file $vcf_File \\
    --output_file VariantAnnotation/VEP/$datasetID/$patientID/$sampleID/vep.vcf.gz \\
    --vcf \\
    --plugin SingleLetterAA \\
    --compress_output gzip \\
    --offline \\
    --cache \\
    --dir_cache $params.vep_dir_cache \\
    --cache_version $params.vep_cache_version \\
    --assembly $params.assembly \\
    --force \\
    --no_stats \\
    --everything \\
    --use_given_ref \\
    --fasta  $params.ref_genome \\
    --fork 23
    
    """
}
