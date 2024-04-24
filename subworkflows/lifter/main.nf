//
// LIFTER SUB-WORKFLOW
//

include { GET_POSITIONS } from "../../modules/mpileup/main_vcf"
include { BCFTOOLS_MPILEUP } from "../../modules/mpileup/main"
include { JOIN_POSITIONS } from "../../modules/mpileup/main_join"


workflow LIFTER {
    take:
        rds
        tumor_bam


    main:

        GET_POSITIONS(rds.groupTuple(by: [0,1]))
        BCFTOOLS_MPILEUP(GET_POSITIONS.out.bed.transpose(by: [2,3]), tumor_bam)
        res = VCF_PROCESSING.out.rds.join(BCFTOOLS_MPILEUP.out.vcf, by: [0,1,2])
        JOIN_POSITIONS(res, GET_POSITIONS.out.pos.transpose(by: 2))

    emit:
        JOIN_POSITIONS.out.rds

}