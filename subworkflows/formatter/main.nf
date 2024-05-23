//
// FORMATTING SUB-WORKFLOW
//

include { CNA_PROCESSING } from '../../modules/cna2CNAqc/main'
include { VCF_PROCESSING } from '../../modules/vcf2CNAqc/main'
include { RDS_PROCESSING } from '../../modules/CNAqc2tsv/main'

workflow FORMATTER {
    take:
        input
        extension
    
    main:

        if (extension == "vcf"){
                out = VCF_PROCESSING(input)
        } else if (extension == "cna"){
                out = CNA_PROCESSING(input)
         } else if (extension == "rds"){
                out = RDS_PROCESSING(input)
         }

    emit:
        out
}
