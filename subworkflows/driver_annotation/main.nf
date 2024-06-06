//
// DRIVER_ANNOTATION SUB-WORKFLOW
//

include { ANNOTATE_DRIVER } from '../../modules/driver_annotation/main'


workflow DRIVER_ANNOTATION {
    take:
        rds
        cancer_type
    
    main:
        ANNOTATE_DRIVER(rds, cancer_type)


    emit:
        ANNOTATE_DRIVER.out.rds

}