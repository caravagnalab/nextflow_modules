// nextflow.enable.moduleBinaries = true

process PYCLONEVI {
    publishDir (
      params.publish_dir,
      mode: "copy",
      pattern: "$patientID/pyclonevi/*"
    )
               
    input:

      tuple val(patientID), path(joint_table)

    output:
      tuple val(patientID), path("$patientID/ctree/ctree_input_pyclonevi.csv"), emit: ctree_input
      tuple val(patientID), path("$patientID/pyclonevi/all_fits.h5"),  emit: pyclone_all_fits
      tuple val(patientID), path("$patientID/pyclonevi/best_fit.tsv"), emit: pyclone_best_fit
      // tuple val(patientID), path("$patientID/pyclonevi/pyclone_report.pdf"), emit: pyclone_report

    script:
      def args = task.ext.args ?: ''
      def n_cluster_arg                    = args.n_cluster                     ?  "$args.n_cluster" : ""
      def density_arg                    = args.density                     ?  "$args.density" : ""
      def n_grid_point_arg                    = args.n_grid_point                     ?  "$args.n_grid_point" : ""
      def n_restarts_arg                    = args.n_restarts                     ?  "$args.n_restarts" : ""


      """
      mkdir -p $patientID/pyclonevi
      mkdir -p $patientID/ctree

      pyclone-vi fit -i $joint_table \\
      -o $patientID/pyclonevi/all_fits.h5 -c $n_cluster_arg \\
      -d $density_arg \\
      --num-grid-points $n_grid_point_arg \\
      --num-restarts $n_restarts_arg

      pyclone-vi write-results-file -i $patientID/pyclonevi/all_fits.h5  \\
      -o $patientID/pyclonevi/best_fit.tsv

      #python3 /orfeo/cephfs/scratch/cdslab/ggandolfi/nextflow_modules/modules/pyclonevi/pyclone_ctree.py $joint_table \\
      python3 $moduleDir/pyclone_ctree.py $joint_table \\
      $patientID/pyclonevi/best_fit.tsv \\
      $patientID/ctree/ctree_input_pyclonevi.csv
      
      #python3 /orfeo/cephfs/scratch/cdslab/ggandolfi/nextflow_modules/modules/pyclonevi/pyclone_plot.py $joint_table \\
      #$patientID/pyclonevi/best_fit.tsv \\
      #$patientID/pyclonevi/pyclone_report.pdf


      """
}
