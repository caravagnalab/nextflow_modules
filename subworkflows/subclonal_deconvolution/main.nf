//
// SUBCLONAL DECONVOLUTION WORKFLOW
//

include { PYCLONEVI } from '../../modules/pyclonevi/main'
include { MOBSTERh } from '../../modules/mobsterh/main'
include { CTREE as CTREE_PYCLONEVI } from '../../modules/ctree/main'
include { CTREE as CTREE_MOBSTERh } from '../../modules/ctree/main'

workflow SUBCLONAL_DECONVOLUTION {
    take: 
    joint_table

    main:
    if (params.tools && params.tools.split(',').contains('pyclone-vi')) {
        PYCLONEVI(joint_table) // run pyclone-vi
        CTREE_PYCLONEVI(PYCLONEVI.out.ctree_input)
        ctree_plot_pyclonevi = CTREE_PYCLONEVI.out.ctree_plot
        emit:
        ctree_plot_pyclonevi
    }

    if (params.tools && params.tools.split(',').contains('mobsterh')) {
      MOBSTERh(joint_table) // run mobsterh
      CTREE_MOBSTERh(MOBSTERh.out.ctree_input)
      ctree_plot_mobsterh = CTREE_MOBSTERh.out.ctree_plot
      emit:
      ctree_plot_mobsterh
    }
}
