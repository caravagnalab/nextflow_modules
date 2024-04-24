//
// FORMATTING SUB-WORKFLOW
//

include { CNA_PROCESSING } from '../../modules/cna2CNAqc/main'
include { VCF_PROCESSING } from '../../modules/vcf2CNAqc/main'


workflow FORMATTER {
    take:
        input
    
    main:

        CNA_PROCESSING(input)
        VCF_PROCESSING(input)

    emit:

}