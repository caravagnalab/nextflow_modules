//
// QC SUB-WORKFLOW
//

//include { TINC } from '../../modules/TINC/main'
include { CNAQC } from '../../modules/CNAqc/main'
include { JOIN_CNAQC } from '../../modules/join_CNAqc/main'


workflow QC {
    take: 
        cna
        vcf

    main:

        //TINC(cna, vcf)
        CNAQC(cna, vcf)
        JOIN_CNAQC(CNAQC.out.qc_rds.groupTuple(by: [0,1]))
    
    emit:
        rds_cnaqc = CNAQC.out.qc_rds
        plot_cnaqc_rds = CNAQC.out.plot_rds
        plot_cnaqc_data = CNAQC.out.plot_pdf_data
        plot_cnaqc_qc = CNAQC.out.plot_pdf_qc

        //plot_rds_tinc = TINC.out.plot_rds
        //rds_tinc = TINC.out.fit_rds
        //pdf_tinc = TINC.out.plot_pdf

        rds_join = JOIN_CNAQC.out.rds
}
