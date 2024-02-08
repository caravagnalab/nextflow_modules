//
// VARIANT ANNOTATION WORKFLOW
//

include { VEP_ANNOTATE } from "${baseDir}/modules/VEP/main"
include { VCF2MAF } from "${baseDir}/modules/vcf2maf/main"
include { MAFTOOLS } from "${baseDir}/modules/MAFTOOLS/main"
include { ANNOVAR_ANNOTATE } from "${baseDir}/modules/annovar/main"
include { ANNOVAR2MAF } from "${baseDir}/modules/annovar2maf/main"

workflow VARIANT_ANNOTATION {
    take: 
    vcf_File

    main:
    if (params.tools && params.tools.split(',').contains('vep')) {
        vep_out = VEP_ANNOTATE(vcf_File) // run VEP
        vcf2maf_out = VCF2MAF(VEP_ANNOTATE.out.vep_output)
        maf_out = MAFTOOLS(VCF2MAF.out.vcf2maf_out.groupTuple(by: 0))
        emit:
        vep_out
        vcf2maf_out
        maf_out
    }

    if (params.tools && params.tools.split(',').contains('annovar')) {
      annovar_out = ANNOVAR_ANNOTATE(vcf_File)
      annovar2maf_out = ANNOVAR2MAF(annovar_out)
      emit:
      annovar_out
      annovar2maf_out
    }
}
