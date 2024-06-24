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
    plot_pdf = null
    plot_rds = null
    signatures_nmfOut = null 
    bestConf = null
    sign_cv = null


    if (params.tools && params.tools.split(',').contains('sparsesignatures')) {
        out = FORMATTER_RDS(joint_table, "rds")
        SPARSE_SIGNATURES(out.groupTuple(by: 0)) // run SparseSignatures
        
        plot_pdf = SPARSE_SIGNATURES.out.signatures_plot_pdf
        plot_rds = SPARSE_SIGNATURES.out.signatures_plot_rds
        signatures_nmfOut = SPARSE_SIGNATURES.out.signatures_nmfOut_rds
        bestConf = SPARSE_SIGNATURES.out.signatures_bestConf_rds
        sign_cv = SPARSE_SIGNATURES.out.signatures_cv_rds
    } 
    emit:
    plot_pdf
    plot_rds
    signatures_nmfOut
    bestConf
    sign_cv

    //if (params.tools && params.tools.split(',').contains('sigprofiler')) {
    //    out = FORMATTER_RDS(joint_table, "rds")
    //    SigProfiler_out = SIGPROFILER(out) // run SigProfiler
    //  
    //    emit:
    //    SigProfiler_out
    //}
}
