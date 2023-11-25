process PYCLONEVI {

    publishDir params.publish_dir, mode: 'copy'

    input:

      tuple val(patientID), path(joint_table)

    output:

      tuple path("$patientID/pyclonevi/all_fits.h5"), path("$patientID/pyclonevi/best_fit.tsv")

    script:
      def args = task.ext.args ?: ''
      def n_cluster_arg                    = args.n_cluster                     ?  "$args.n_cluster" : ""
      def density_arg                    = args.density                     ?  "$args.density" : ""
      def n_grid_point_arg                    = args.n_grid_point                     ?  "$args.n_grid_point" : ""
      def n_restarts_arg                    = args.n_restarts                     ?  "$args.n_restarts" : ""


      """
      mkdir -p $patientID/pyclonevi

      pyclone-vi fit -i $joint_table \\
      -o $patientID/pyclonevi/all_fits.h5 -c $n_cluster_arg \\
      -d $density_arg \\
      --num-grid-points $n_grid_point_arg \\
      --num-restarts $n_restarts_arg

      pyclone-vi write-results-file -i $patientID/pyclonevi/all_fits.h5  \\
      -o $patientID/pyclonevi/best_fit.tsv

      """
}
