//
// CNAqc WORKFLOW
//

include { CNA_PROCESSING } from '../../modules/cna2CNAqc/main'
include { VCF_PROCESSING } from '../../modules/vcf2CNAqc/main'
include { CNAQC } from '../../modules/CNAqc/main'
include { JOIN_CNAQC } from '../../modules/join_CNAqc/main'

include { GET_POSITIONS } from "../../modules/mpileup/main_vcf"
include { BCFTOOLS_MPILEUP } from "../../modules/mpileup/main"
include { JOIN_POSITIONS } from "../../modules/mpileup/main_join"


workflow QC {
    take: 
    vcf
    cnv
    tumor_bam

    main:

    CNA_PROCESSING(cnv)
    VCF_PROCESSING(vcf)

    GET_POSITIONS(VCF_PROCESSING.out.rds.groupTuple(by: [0,1]))
    BCFTOOLS_MPILEUP(GET_POSITIONS.out.bed.transpose(by: [2,3]), tumor_bam)


    ch = VCF_PROCESSING.out.rds.join(BCFTOOLS_MPILEUP.out.vcf, by: [0,1,2])
    JOIN_POSITIONS(ch, GET_POSITIONS.out.pos.transpose(by: 2))

    //CNAQC(CNA_PROCESSING.out.rds, JOIN_POSITIONS.out.rds)
    //JOIN_CNAQC(CNAQC.out.rds.groupTuple(by: [0,1]))
    
    emit:
    //BCFTOOLS_MPILEUP.out
    //GET_POSITIONS.out.pos
    JOIN_POSITIONS.out

    //CNAQC.out.rds
    //CNAQC.out.pdf

    //JOIN_CNAQC.out.rds
}
