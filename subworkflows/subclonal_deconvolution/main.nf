//
// SUBCLONAL DECONVOLUTION WORKFLOW
//
include { VIBER as VIBER_SINGLE } from '../../modules/viber/main'
include { VIBER as VIBER_MULTI } from '../../modules/viber/main'
include { PYCLONEVI as PYCLONEVI_MULTI} from '../../modules/pyclonevi/main'
include { PYCLONEVI as PYCLONEVI_SINGLE} from '../../modules/pyclonevi/main'
include { CTREE as CTREE_PYCLONEVI } from '../../modules/ctree/main'
include { MOBSTERh as MOBSTERh_SINGLE } from '../../modules/mobsterh/main'
include { MOBSTERh as MOBSTERh_MULTI } from '../../modules/mobsterh/main'
include { JOINT_FIT } from '../../modules/joint_fit/main'
include { CTREE as CTREE_MOBSTERh } from '../../modules/ctree/main'

workflow SUBCLONAL_DECONVOLUTION {
    take: 
    joint_table // grouped by patient

    main:
    // single sample subclonal deconvolution

    if (params.step && params.step.split(',').contains('subclonal_singlesample')){
        input_joint_table=joint_table.transpose(by: [2]) // split by patient to have list of samples
        // viber single sample
        if (params.tools && params.tools.split(',').contains('viber')){ 
            t=VIBER_SINGLE(input_joint_table)
            //CTREE_VIBER(VIBER_SINGLE.out.ctree_input)
            //ctree_plot_viber = CTREE_VIBER.out.ctree_plot
            emit:
            t
            //ctree_plot_viber
        // pyclone-vi single sample
        }
	if (params.tools && params.tools.split(',').contains('pyclone-vi')){
            // need to parse m_cnaqc obj to convert to tsv
            //CNAQC_TO_TSV(input_joint_table) // list of patient and sample ids
            t=PYCLONEVI_SINGLE(input_joint_table) 
            //CTREE_PYCLONEVI(PYCLONEVI_SINGLE.out.ctree_input)
            //ctree_plot_pyclonevi = CTREE_PYCLONEVI.out.ctree_plot
            emit:
            t
            //ctree_plot_pyclonevi
        }
        if (params.tools && params.tools.split(',').contains('mobster')){
            t=MOBSTERh_SINGLE(input_joint_table)
            //CTREE_PYCLONEVI(PYCLONEVI_SINGLE.out.ctree_input)
            //ctree_plot_pyclonevi = CTREE_PYCLONEVI.out.ctree_plot
            emit:
            t
            //ctree_plot_pyclonevi
        }
    }
    //multi sample subclonal deconvolution

    if (params.step && params.step.split(',').contains('subclonal_multisample')){
        if (params.tools && params.tools.split(',').contains('mobster')){
            input_joint_table=joint_table.transpose(by: [1])
            MOBSTERh_MULTI(input_joint_table) // tuple val(patientID), val(sampleID), path("$outDir/mobster_joint*"), emit: mobster_joint
            input_joint_table=JOINT_FIT((MOBSTERh_MULTI.out.mobster_joint).groupTuple(by: [0]))
            }
        if (params.tools && params.tools.split(',').contains('viber')){
            t=VIBER_MULTI(joint_table)
            //CTREE_VIBER(VIBER_MULTI.out.ctree_input)
            //ctree_plot_viber = CTREE_VIBER.out.ctree_plot
            emit:
            t
            //ctree_plot_viber
            }
        if (params.tools && params.tools.split(',').contains('pyclone-vi')){
            CNAQC_TO_TSV(input_joint_table) // list of patient and sample ids
            t=PYCLONEVI_MULTI(CNAQC_TO_TSV.out.input_table)
            //CTREE_PYCLONEVI(PYCLONEVI_MULTI.out.ctree_input)
            //ctree_plot_pyclonevi = CTREE_PYCLONEVI.out.ctree_plot
            emit:
            t
            //ctree_plot_pyclonevi
            }
    }
}
