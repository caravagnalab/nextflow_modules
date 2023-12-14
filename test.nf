#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process pullDockerImage {
    container 'docker'
    
    script:
    """
    singularity pull docker://bdgenomics/rhapsody
    """
}

workflow {
    pullDockerImage()
}
