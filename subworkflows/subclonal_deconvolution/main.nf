//
// SUBCLONAL DECONVOLUTION WORKFLOW
//

include { PYCLONEVI } from '../../modules/pyclonevi/main'
include { CTREE } from '../../modules/ctree/main'

workflow SUBCLONAL_DECONVOLUTION {
    take: 
    joint_table // input for pyclone

    main:
    if (params.tools && params.tools.split(',').contains('pyclone-vi')) {
    PYCLONEVI(joint_table) // run pyclone-vi
    CTREE(PYCLONEVI.out.ctree_input)
    ctree_plot = CTREE.out.ctree_plot
    emit:

    ctree_plot
    }
}
