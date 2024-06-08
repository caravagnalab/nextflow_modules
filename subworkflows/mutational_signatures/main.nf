//
// MUTATIONAL SIGNATURES DECONVOLUTION WORKFLOW
//

include { FORMATTER as FORMATTER_RDS} from "${baseDir}/nextflow_modules/subworkflows/formatter/main"
include { SPARSE_SIGNATURES } from "${baseDir}/nextflow_modules/modules/SparseSignatures/main"
//include { SIGPROFILER } from "${baseDir}/nextflow_modules/modules/SigProfiler/main"


workflow MUTATIONAL_SIGNATURES {
    take: 
    join_cnaqc

    main:
   
    
    if (params.tools && params.tools.split(',').contains('sparsesignatures')) {
        out = FORMATTER_RDS(join_cnaqc, "rds")
        SparseSig_out = SPARSE_SIGNATURES(out) // run SparseSignatures
        
        emit:
        SparseSig_out
    }

    if (params.tools && params.tools.split(',').contains('sigprofiler')) {
        out = FORMATTER_RDS(join_cnaqc, "rds")
        SigProfiler_out = SIGPROFILER(out) // run SigProfiler
      
        emit:
        SigProfiler_out
    }
}
