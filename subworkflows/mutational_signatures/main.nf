//
// MUTATIONAL SIGNATURES WORKFLOW
//

include { SIG_PROFILER } from "${baseDir}/../modules/SigProfiler/main"
include { SPARSE_SIG } from "${baseDir}/../modules/SparseSig/main"


workflow MUTATIONAL_SIGNATURES {
    take:
    tsv_joint

    main:
    if (params.tools && params.tools.split(',').contains('SigProfiler')) {
        SigProfiler_out = SIG_PROFILER(tsv_joint) // run SigProfiler
   
        emit:
        SigProfiler_out
        
       
    }

    if (params.tools && params.tools.split(',').contains('SparseSig')) {
      SparseSig_out = SPARSE_SIGNATURES(tsv_joint)
  
      emit:
      SparseSig_out
    }
}
