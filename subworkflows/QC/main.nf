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

    // CNA_PROCESSING(cnv)
    // VCF_PROCESSING(vcf)

    // GET_POSITIONS(VCF_PROCESSING.out.rds.groupTuple(by: [0,1]))
    // BCFTOOLS_MPILEUP(GET_POSITIONS.out.bed.transpose(by: [2,3]), tumor_bam)

    // ch = VCF_PROCESSING.out.rds.join(BCFTOOLS_MPILEUP.out.vcf, by: [0,1,2])
    // JOIN_POSITIONS(ch, GET_POSITIONS.out.pos.transpose(by: 2))

    CNAQC(CNA_PROCESSING.out.rds, JOIN_POSITIONS.out.rds)
    JOIN_CNAQC(CNAQC.out.rds.groupTuple(by: [0,1]))
    
    emit:
    //CNAQC.out.rds
    //CNAQC.out.pdf

    //JOIN_CNAQC.out.rds
}
