//
// VARIANT ANNOTATION WORKFLOW
//

include { VEP_ANNOTATE } from "${baseDir}/modules/VEP/main"
include { VCF2MAF } from "${baseDir}/modules/vcf2maf/main"
include { MAFTOOLS } from "${baseDir}/modules/MAFTOOLS/main"


workflow VARIANT_ANNOTATION {
    take: 
        vcf

    main:
        VEP_ANNOTATE(vcf)
        VCF2MAF(VEP_ANNOTATE.out.vep_output)
        MAFTOOLS(VCF2MAF.out.vcf2maf_out.groupTuple(by: 0))

    emit:
        vep = VEP_ANNOTATE.out.vep_output 
        vcf2maf = VCF2MAF.out.vcf2maf_out
        
        maf_rds = MAFTOOLS.out.maf_results
        maf_oncoplot = MAFTOOLS.out.oncoplot
        maf_summary = MAFTOOLS.out.summary_plot

}
