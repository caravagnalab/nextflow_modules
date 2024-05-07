//
// QC SUB-WORKFLOW
//

include { CNAQC } from '../../modules/CNAqc/main'
include { JOIN_CNAQC } from '../../modules/join_CNAqc/main'


workflow QC {
    take: 
        cna
        vcf

    main:

        CNAQC(CNA_PROCESSING.out.rds, JOIN_POSITIONS.out.rds)
        JOIN_CNAQC(CNAQC.out.rds.groupTuple(by: [0,1]))
    
    emit:
        CNAQC.out.rds
        CNAQC.out.pdf

        JOIN_CNAQC.out.rds
}
