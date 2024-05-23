//
// MUTATIONAL SIGNATURES DECONVOLUTION WORKFLOW
//

include { RDS_PROCESSING } from "${baseDir}/modules/CNAqc2tsv/main"
include { SPARSE_SIGNATURES } from "${baseDir}/modules/SparseSignatures/main"
include { SIGPROFILER } from "${baseDir}/modules/SigProfiler/main"


workflow MUTATIONAL_SIGNATURES {
    take: 
    joint_table

    main:
    
    { rds_processing_out = RDS_PROCESSING(join_cnaqc) 

        emit:
        rds_processing_out
    }

    
    
    if (params.tools && params.tools.split(',').contains('SparseSignatures')) {
        SparseSig_out = SPARSESIGNATURES(rds_processing_out) // run SparseSignatures
        
        emit:
        SparseSig_out
    }

    if (params.tools && params.tools.split(',').contains('SigProfiler')) {
      SigProfiler_out = SIGPROFILER(rds_processing_out) // run SigProfiler
      
      emit:
      SigProfiler_out
    }
}
