process PLATYPUS_CALL_VARIANTS {

    publishDir params.publish_dir, mode: 'copy'

    input:

    tuple val(patientID), path(tumor_bamFile), path(tumor_baiFile), path(normal_bamFile), path(normal_baiFile), path(snv_vcfFile), path(snv_tbiFile), path(indel_vcfFile), path(indel_tbiFile) 

    output:

    path("$patientID/Platypus/*.vcf")

    script:

    """
    #!/bin/bash

    mkdir -p $patientID/Platypus

    vcfs=\$(echo $snv_vcfFile,$indel_vcfFile | tr ' ' ',')
    bams=\$(echo $tumor_bamFile,$normal_bamFile | tr ' ' ',')
    
    platypus callVariants \\
    --refFile=$params.ref_genome \\
    --bamFiles=\$bams \\
    --output=$patientID/Platypus/Platypus_${patientID}_joint.vcf \\
    --source=\$vcfs \\
    --filterReadPairsWithSmallInserts=0 \\
    --maxReads=100000000 \\
    --maxVariants=100 \\
    --minPosterior=0 \\
    --nCPU=23 \\
    --getVariantsFromBAMs=0
    """
}
