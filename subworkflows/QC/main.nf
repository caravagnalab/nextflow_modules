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

        CNAQC(cna, vcf)
        JOIN_CNAQC(CNAQC.out.rds.groupTuple(by: [0,1]))
    
    emit:
        rds_cnaqc = CNAQC.out.rds
        pdf_cnaqc = CNAQC.out.pdf

        rds_join = JOIN_CNAQC.out.rds
}
