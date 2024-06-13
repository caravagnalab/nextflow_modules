//
// MUTATIONAL SIGNATURES DECONVOLUTION WORKFLOW
//

include { FORMATTER as FORMATTER_RDS} from "../../subworkflows/formatter/main"
include { SPARSE_SIGNATURES } from "../../modules/SparseSignatures/main"
include { SIG_PROFILER } from "../../modules/SigProfiler/main"


workflow MUTATIONAL_SIGNATURES {
    take: 
    joint_table

    main:
    
    if (params.tools && params.tools.split(',').contains('sparsesignatures')) {
        out = FORMATTER_RDS(joint_table, "rds")
        SparseSig_out = SPARSE_SIGNATURES(out) // run SparseSignatures
        
        emit:
        SparseSig_out
    }

    if (params.tools && params.tools.split(',').contains('sigprofiler')) {
        out = FORMATTER_RDS(joint_table, "rds")
        SigProfiler_out = SIGPROFILER(out) // run SigProfiler
      
        emit:
        SigProfiler_out
    }
}
