//
// Function to initialise default values and to generate a Groovy Map of available options for nf-core modules
//
def initOptions(Map args) {
    def Map options = [:]
    options.args = args.args ?: ''
    options.norm_method = args.norm_method ?: "median"  // normalization.method
    options.window = args.window ?: 1e5  // window
    options.gamma = args.gamma ?: 280  // gamma
    options.kmin = args.kmin ?: 300  // kmin
    options.min_reads_baf = args.min_reads_baf ?: 50
    options.min_reads = args.min_reads ?: 50
    options.min_reads_normal = args.min_reads_normal ?: 15
    options.max_mut_types = args.max_mut_types ?: 1

    options.low_cell = args.low_cell ?: 1
    options.up_cell = args.up_cell ?: 1
    options.low_ploidy = args.low_ploidy ?: 1
    options.up_ploidy = args.up_ploidy ?: 1
    options.is_female = args.is_female ?: "False"

    options.sample_id = args.sample_id ?: ""
    options.out_dir = args.out_dir ?: ""
    // options.publish_by_meta = args.publish_by_meta ?: []
    // options.publish_dir     = args.publish_dir ?: ''
    // options.publish_files   = args.publish_files
    // options.suffix          = args.suffix ?: ''
    return options
}
