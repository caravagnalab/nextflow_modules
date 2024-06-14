//
// MUTATIONAL SIGNATURES DECONVOLUTION WORKFLOW
//

include { FORMATTER as FORMATTER_RDS} from "../../subworkflows/formatter/main"
include { SPARSE_SIGNATURES } from "../../modules/SparseSignatures/main"
//include { SIG_PROFILER } from "../../modules/SigProfiler/main"


workflow MUTATIONAL_SIGNATURES {
    take: 
    joint_table

    main:
    
    //if (params.tools && params.tools.split(',').contains('sparsesignatures')) {
    out = FORMATTER_RDS(joint_table, "rds")
    SPARSE_SIGNATURES(out) // run SparseSignatures
    
    emit:
    plot_pdf = SPARSE_SIGNATURES.out.signatures_plot_pdf
    signatures_nmfOut = SPARSE_SIGNATURES.out.signatures_nmfOut_rds
    plot_rds = SPARSE_SIGNATURES.out.signatures_plot_rds
    bestConf = SPARSE_SIGNATURES.out.signatures_bestConf_rds
    sign = SPARSE_SIGNATURES.out.signatures_cv_rds
    //}

    //if (params.tools && params.tools.split(',').contains('sigprofiler')) {
    //    out = FORMATTER_RDS(joint_table, "rds")
    //    SigProfiler_out = SIGPROFILER(out) // run SigProfiler
    //  
    //    emit:
    //    SigProfiler_out
    //}
}
