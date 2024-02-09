//
// VARIANT ANNOTATION WORKFLOW
//

include { VEP_ANNOTATE } from "${baseDir}/modules/VEP/main"
include { VCF2MAF } from "${baseDir}/modules/vcf2maf/main"
include { MAFTOOLS } from "${baseDir}/modules/MAFTOOLS/main"


workflow VARIANT_ANNOTATION {
    take: 
    vcf_File

    main:
        vep_out = VEP_ANNOTATE(vcf_File)
        vcf2maf_out = VCF2MAF(VEP_ANNOTATE.out.vep_output)
        maf_out = MAFTOOLS(VCF2MAF.out.vcf2maf_out.groupTuple(by: 0))

        emit:
        vep_out
        vcf2maf_out
        maf_out

}
