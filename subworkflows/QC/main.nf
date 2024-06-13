//
// QC SUB-WORKFLOW
//

include { TINC } from '../../modules/TINC/main'
include { CNAQC } from '../../modules/CNAqc/main'
include { JOIN_CNAQC } from '../../modules/join_CNAqc/main'


workflow QC {
    take: 
        cna
        vcf

    main:

        TINC(cna, vcf)
        CNAQC(cna, vcf)
        JOIN_CNAQC(CNAQC.out.rds.groupTuple(by: [0,1]))
    
    emit:

        plot_rds_tinc = TINC.out.plot_rds
        rds_tinc = TINC.out.fit_rds
        pdf_tinc = TINC.out.plot_pdf

        rds_cnaqc = CNAQC.out.rds
        pdf_cnaqc = CNAQC.out.pdf

        rds_join = JOIN_CNAQC.out.rds
}
