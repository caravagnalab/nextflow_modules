//
// SUBCLONAL DECONVOLUTION WORKFLOW
//
include { VIBER as VIBER_SINGLE } from "../../modules/viber/main"
include { VIBER as VIBER_MULTI } from "../../modules/viber/main"
include { PYCLONEVI as PYCLONEVI_MULTI} from "../../modules/pyclonevi/main"
include { PYCLONEVI as PYCLONEVI_SINGLE} from "../../modules/pyclonevi/main"
include { CTREE as CTREE_PYCLONEVI } from "../../modules/ctree/main"
include { CTREE as CTREE_VIBER } from "../../modules/ctree/main"
include { CTREE as CTREE_MOBSTERh } from "../../modules/ctree/main"
include { MOBSTERh as MOBSTERh_SINGLE } from "../../modules/mobsterh/main"
include { MOBSTERh as MOBSTERh_MULTI } from "../../modules/mobsterh/main"
include { FORMATTER as FORMATTER_RDS_SINGLE} from "../../subworkflows/formatter/main"
include { FORMATTER as FORMATTER_RDS_MULTI} from "../../subworkflows/formatter/main"
include { JOINT_FIT } from "../../modules/joint_fit/main"

workflow SUBCLONAL_DECONVOLUTION {
    take: 
    joint_table // grouped by patient

    main:
    
    // single sample subclonal deconvolution
    if (params.mode && params.mode.split(",").contains("singlesample")) {
        input_joint_table = joint_table.transpose(by: [2]) // split by patient to have list of samples

        // viber single sample
        if (params.tools && params.tools.split(",").contains("viber")) { 
            t = VIBER_SINGLE(input_joint_table)
            // working but not producing outputs with the current VIBER function `get_clone_trees()` when only one sample
            CTREE_VIBER(VIBER_SINGLE.out.viber_rds)
            emit:
            t
        }
        if (params.tools && params.tools.split(",").contains("pyclone-vi")) {
            // need to parse m_cnaqc obj to convert to tsv
            // CNAQC_TO_TSV(input_joint_table) // list of patient and sample ids
            FORMATTER_RDS_SINGLE(input_joint_table, "rds")
            t = PYCLONEVI_SINGLE(FORMATTER_RDS_SINGLE.out) 
            // CTREE_PYCLONEVI(PYCLONEVI_SINGLE.out.ctree_input)
            emit:
            t
        }
        if (params.tools && params.tools.split(",").contains("mobster")) {
            t = MOBSTERh_SINGLE(input_joint_table)
            CTREE_MOBSTERh(MOBSTERh_SINGLE.out.mobster_best_rds)
            emit:
            t
        }
    }


    // multi sample subclonal deconvolution
    if (params.mode && params.mode.split(",").contains("multisample")) {

        if (params.tools && params.tools.split(",").contains("mobster")) {
            // run MOBSTER on all samples independently
            joint_table_singlesample = joint_table.transpose(by: [2])
            // tuple val(datasetID), val(patientID), val(sampleID), path("${outDir}/*/mobsterh_st_best_fit.rds"), emit: mobster_best_rds
            MOBSTERh_MULTI(joint_table_singlesample)
            // collect all results and group by patient
            joint_table = JOINT_FIT( (MOBSTERh_MULTI.out.mobster_best_rds).groupTuple(by: [1]) )
        }

        if (params.tools && params.tools.split(",").contains("viber")) {
            t = VIBER_MULTI(joint_table)
            CTREE_VIBER(VIBER_MULTI.out.viber_rds)
            emit:
            t
        }
        if (params.tools && params.tools.split(",").contains("pyclone-vi")) {
            // CNAQC_TO_TSV(joint_table) // list of patient and sample ids
            FORMATTER_RDS_MULTI(joint_table, "rds")
            t = PYCLONEVI_MULTI(FORMATTER_RDS_MULTI.out)
            // CTREE_PYCLONEVI(PYCLONEVI_MULTI.out.ctree_input)
            emit:
            t
        }
    }
}
