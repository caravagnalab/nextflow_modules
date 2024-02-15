//
// CNAqc WORKFLOW
//

include { CNA_PROCESSING } from '../../modules/cna2CNAqc/main'
include { VCF_PROCESSING } from '../../modules/vcf2CNAqc/main'
include { CNAQC_ANALYSIS } from '../../modules/CNAqc/main'
include { JOIN_CNAQC } from '../../modules/join_CNAqc/main'

workflow CNAQC {
    take: 
    vcf
    input_CNA

    main:

    CNA_PROCESSING(input_CNA) 
    VCF_PROCESSING(vcf)
    CNAQC_ANALYSIS(CNA_PROCESSING.out.rds, VCF_PROCESSING.out.rds)

    JOIN_CNAQC(CNAQC_ANALYSIS.out.rds.groupTuple(by: [0,1]))
    
    emit:

    CNAQC_ANALYSIS.out.rds
    CNAQC_ANALYSIS.out.pdf

    JOIN_CNAQC.out.rds
    JOIN_CNAQC.out.mut_tsv
    JOIN_CNAQC.out.cna_tsv
}
