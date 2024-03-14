//
// CNAqc WORKFLOW
//

include { CNA_PROCESSING } from '../../modules/cna2CNAqc/main'
include { VCF_PROCESSING } from '../../modules/vcf2CNAqc/main'
include { CNAQC } from '../../modules/CNAqc/main'
include { JOIN_CNAQC } from '../../modules/join_CNAqc/main'

workflow QC {
    take: 
    vcf
    cnv

    main:

    CNA_PROCESSING(cnv)
    VCF_PROCESSING(vcf)
    CNAQC(CNA_PROCESSING.out.rds, VCF_PROCESSING.out.rds)

    JOIN_CNAQC(CNAQC.out.rds.groupTuple(by: [0,1]))
    
    emit:

    CNAQC.out.rds
    CNAQC.out.pdf

    JOIN_CNAQC.out.rds
}
