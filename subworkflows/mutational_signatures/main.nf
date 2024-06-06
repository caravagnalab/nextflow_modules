//
// MUTATIONAL SIGNATURES DECONVOLUTION WORKFLOW
//

include { SPARSE_SIGNATURES } from "${baseDir}/modules/SparseSignatures/main"
include { SIGPROFILER } from "${baseDir}/modules/SigProfiler/main"


workflow MUTATIONAL_SIGNATURES {
    take: 
    joint_table

    main:
    if (params.tools && params.tools.split(',').contains('SparseSignatures')) {
        SparseSig_out = SPARSESIGNATURES(joint_table) // run SparseSignatures
        
        emit:
        SparseSig_out
    }

    if (params.tools && params.tools.split(',').contains('SigProfiler')) {
      SigProfiler_out = SIGPROFILER(joint_table) // run SigProfiler
      
      emit:
      SigProfiler_out
    }
}
