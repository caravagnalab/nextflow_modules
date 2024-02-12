//
// SUBCLONAL DECONVOLUTION WORKFLOW
//
include { VIBER as VIBER_SINGLE } from '../../modules/viber/main'
include { VIBER as VIBER_MULTI } from '../../modules/viber/main'
include { PYCLONEVI as PYCLONEVI_MULTI} from '../../modules/pyclonevi/main'
include { PYCLONEVI as PYCLONEVI_SINGLE} from '../../modules/pyclonevi/main'
include { CTREE as CTREE_PYCLONEVI } from '../../modules/ctree/main'
//include { MOBSTERh } from '../../modules/mobsterh/main'
//include { CTREE as CTREE_MOBSTERh } from '../../modules/ctree/main'

workflow SUBCLONAL_DECONVOLUTION {
    take: 
    joint_table

    main:
    // single sample subclonal deconvolution

    if (params.step && params.step.split(',').contains('subclonal_singlesample')){
        input_joint_table=joint_table.transpose(by: [1]) 
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
            t=PYCLONEVI_SINGLE(input_joint_table)
            //CTREE_PYCLONEVI(PYCLONEVI_SINGLE.out.ctree_input)
            //ctree_plot_pyclonevi = CTREE_PYCLONEVI.out.ctree_plot
            emit:
            t
            //ctree_plot_pyclonevi
        }
    }

    //multi sample subclonal deconvolution

    if (params.step && params.step.split(',').contains('subclonal_multisample')){
        //viber multi sample
        if (params.tools && params.tools.split(',').contains('viber')){
            t=VIBER_MULTI(joint_table)
            //CTREE_VIBER(VIBER_MULTI.out.ctree_input)
            //ctree_plot_viber = CTREE_VIBER.out.ctree_plot
            emit:
            t
            //ctree_plot_viber
        // pyclone-vi multi sample
        }
	if (params.tools && params.tools.split(',').contains('pyclone-vi')){
            t=PYCLONEVI_MULTI(joint_table)
            //CTREE_PYCLONEVI(PYCLONEVI_MULTI.out.ctree_input)
            //ctree_plot_pyclonevi = CTREE_PYCLONEVI.out.ctree_plot
            emit:
            t
            //ctree_plot_pyclonevi
        }
    } 
}
