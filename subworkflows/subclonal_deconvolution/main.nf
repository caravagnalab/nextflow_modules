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
    pyclone_fits = null 
    pyclone_best = null
    ctree_pyclone_pdf = null
    viber_pdf = null
    ctree_viber_pdf = null
    mobster_pdf = null
    ctree_mobster_pdf = null

    // single sample subclonal deconvolution
    if (params.mode && params.mode.split(",").contains("singlesample")) {
        input_joint_table = joint_table.transpose(by: [2]) // split by patient to have list of samples

        // viber single sample
        if (params.tools && params.tools.split(",").contains("viber")) { 
            VIBER_SINGLE(input_joint_table)
            // working but not producing outputs with the current VIBER function `get_clone_trees()` when only one sample
            CTREE_VIBER(VIBER_SINGLE.out.viber_rds)
            //emit:
            viber_pdf = VIBER_SINGLE.out.viber_report_pdf
            ctree_viber_pdf = CTREE_VIBER.out.ctree_report_pdf
        }


        if (params.tools && params.tools.split(",").contains("pyclone-vi")) {
            // need to parse m_cnaqc obj to convert to tsv
            // CNAQC_TO_TSV(input_joint_table) // list of patient and sample ids
            FORMATTER_RDS_SINGLE(input_joint_table, "rds")
            PYCLONEVI_SINGLE(FORMATTER_RDS_SINGLE.out) 
            CTREE_PYCLONEVI(PYCLONEVI_SINGLE.out.ctree_input)
            //emit:
            pyclone_fits = PYCLONEVI_SINGLE.out.pyclone_all_fits
            pyclone_best = PYCLONEVI_SINGLE.out.pyclone_best_fit
            ctree_pyclone_pdf = CTREE_PYCLONEVI.out.ctree_report_pdf
        }


        if (params.tools && params.tools.split(",").contains("mobster")) {
            MOBSTERh_SINGLE(input_joint_table)
            CTREE_MOBSTERh(MOBSTERh_SINGLE.out.mobster_best_rds)
            //emit:
            mobster_pdf = MOBSTERh_SINGLE.out.mobster_report_pdf
            ctree_mobster_pdf = CTREE_MOBSTERh.out.ctree_report_pdf
        }
    } 


    // multi sample subclonal deconvolution
    if (params.mode && params.mode.split(",").contains("multisample")) {

        if (params.tools && params.tools.split(",").contains("mobster")) {
            // run MOBSTER on all samples independently
            joint_table_singlesample = joint_table.transpose(by: [2])
            // tuple val(datasetID), val(patientID), val(sampleID), path("${outDir}/*/mobsterh_st_best_fit.rds"), emit: mobster_best_rds
            MOBSTERh_MULTI(joint_table_singlesample)

            input_joint_fit = joint_table_singlesample.join(MOBSTERh_MULTI.out.mobster_best_rds, by: [0,1,2]).groupTuple(by: [0,1,3]) // group by datasetID, patientID, joint_table
            // collect all results and group by patient
            input_joint_table = JOINT_FIT(input_joint_fit)

        } else {
            input_joint_table = joint_table
        }

        if (params.tools && params.tools.split(",").contains("viber")) {
            VIBER_MULTI(input_joint_table)
            CTREE_VIBER(VIBER_MULTI.out.viber_rds)
            //emit:
            viber_pdf = VIBER_MULTI.out.viber_report_pdf
            ctree_viber_pdf = CTREE_VIBER.out.ctree_report_pdf
        }


        if (params.tools && params.tools.split(",").contains("pyclone-vi")) {
            // CNAQC_TO_TSV(joint_table) // list of patient and sample ids
            FORMATTER_RDS_MULTI(input_joint_table, "rds")

            PYCLONEVI_MULTI(FORMATTER_RDS_MULTI.out)
            CTREE_PYCLONEVI(PYCLONEVI_MULTI.out.ctree_input)
            //emit:
            pyclone_fits = PYCLONEVI_MULTI.out.pyclone_all_fits
            pyclone_best = PYCLONEVI_MULTI.out.pyclone_best_fit
            ctree_pyclone_pdf = CTREE_PYCLONEVI.out.ctree_report_pdf
        }
    }

    emit:
    pyclone_fits
    pyclone_best
    ctree_pyclone_pdf
    viber_pdf
    ctree_viber_pdf
    mobster_pdf
    ctree_mobster_pdf
}
