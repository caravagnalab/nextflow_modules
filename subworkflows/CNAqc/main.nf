//
// CNAqc WORKFLOW
//

include { CNA_PROCESSING } from '../../modules/cna2CNAqc/main'
include { VCF_PROCESSING } from '../../modules/vcf2CNAqc/main'
include { CNAQC_ANALYSIS } from '../../modules/CNAqc/main'

workflow CNAQC {
    take: 
    vep_output 
    input_CNA

    main:
    
    CNA_PROCESSING(input_CNA) 
    VCF_PROCESSING(vep_output)
    CNAQC_ANALYSIS(CNA_PROCESSING.out.rds, VCF_PROCESSING.out.rds)
    
    emit:

    CNAQC_ANALYSIS.out.rds
    CNAQC_ANALYSIS.out.pdf
    
}