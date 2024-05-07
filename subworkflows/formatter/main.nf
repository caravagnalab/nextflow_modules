//
// FORMATTING SUB-WORKFLOW
//

include { CNA_PROCESSING } from '../../modules/cna2CNAqc/main'
include { VCF_PROCESSING } from '../../modules/vcf2CNAqc/main'


workflow FORMATTER {
    take:
        input
    
    main:

        VCF_PROCESSING(result.vcf)
        CNA_PROCESSING(result.cna)

    emit:
    
        CNA_PROCESSING.out.rds
        VCF_PROCESSING.out.rds

}