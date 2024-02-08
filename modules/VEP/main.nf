//
// Variant Annotation using VEP
//

process VEP_ANNOTATE {
    publishDir params.publish_dir, mode: 'copy'

    input:

      tuple val(datasetID), val(patientID), val(sampleID), path(vcf_File)  

    output:

      tuple val(datasetID), val(patientID), val(sampleID), path("$datasetID/$patientID/$sampleID/VEP/*.vcf.gz"), emit: vep_output 

    script:

    """
 // Create output directory

    mkdir -p $datasetID/$patientID/$sampleID/VEP

 // Run VEP

    vep --input_file $vcf_File \\
    --output_file $datasetID/$patientID/$sampleID/VEP/vep.vcf.gz \\
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
