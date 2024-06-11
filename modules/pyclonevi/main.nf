process PYCLONEVI {
    publishDir (
      params.publish_dir,
      mode: "copy"
    )
               
    input:

      tuple val(datasetID), val(patientID), val(sampleID), path(joint_table) // from the formatter output

    output:
      tuple val(datasetID), val(patientID), val(sampleID), path("${outDir_ctree}/ctree_input_pyclonevi.csv"), emit: ctree_input
      tuple val(datasetID), val(patientID), val(sampleID), path("${outDir}/all_fits.h5"), emit: pyclone_all_fits
      tuple val(datasetID), val(patientID), val(sampleID), path("${outDir}/best_fit.txt"), emit: pyclone_best_fit
    
    script:
      def args = task.ext.args ?: ''
      def n_cluster_arg                    = args.n_cluster                     ?  "$args.n_cluster" : ""
      def density_arg                    = args.density                     ?  "$args.density" : ""
      def n_grid_point_arg                    = args.n_grid_point                     ?  "$args.n_grid_point" : ""
      def n_restarts_arg                    = args.n_restarts                     ?  "$args.n_restarts" : ""
      def mode                    = args.mode                     ?  "$args.mode" : ""

      if (mode == "singlesample") {
        sampleID_string = sampleID
        outDir = "subclonal_deconvolution/pyclonevi/$datasetID/$patientID/$sampleID"
        outDir_ctree = "subclonal_deconvolution/ctree/$datasetID/$patientID/$sampleID"
      } else if (mode == "multisample"){
        if (!(sampleID instanceof String)) {
          sampleID_string = sampleID.join(" ")
        } else {
          sampleID_string = sampleID
        }
        outDir = "subclonal_deconvolution/pyclonevi/$datasetID/$patientID"
        outDir_ctree = "subclonal_deconvolution/ctree/$datasetID/$patientID"
      }

      all_fits = "${outDir}/all_fits.h5"
      best_fit = "${outDir}/best_fit.txt"
      path_ctree = "${outDir_ctree}/ctree_input_pyclonevi.csv"
      pyclone_joint = "${outDir}/pyclone_joint.tsv"

      """

      mkdir -p $outDir
      mkdir -p $outDir_ctree
      
      # format the input table in order to be pyclone compliant
      python3 $moduleDir/pyclone_utils.py create_pyclone_input $joint_table $patientID $outDir/pyclone_input_all_samples.tsv

      colnames=\$(head -n 1 $outDir/pyclone_input_all_samples.tsv)

      column_number=\$(echo -e "\$colnames" | awk -v col_name="sample_id" 'BEGIN { FS = "\t" } {
          for (i = 1; i <= NF; i++) {
              if (\$i == col_name) {
                  print i
                  exit
              }
          }
          exit 1  # Exit with error if column not found
      }')

      echo -e "\$colnames" > $outDir/pyclone_input.tsv
      for i in $sampleID_string;
        do awk '\$'\$column_number' == "'"\$i"'"' $outDir/pyclone_input_all_samples.tsv >> $outDir/pyclone_input.tsv;
      done
      
      pyclone-vi fit -i $outDir/pyclone_input.tsv -o $all_fits -c $n_cluster_arg -d $density_arg --num-grid-points $n_grid_point_arg --num-restarts $n_restarts_arg
      pyclone-vi write-results-file -i $all_fits -o $best_fit
      python3 $moduleDir/pyclone_utils.py update_hd5_file_with_csv $all_fits $outDir/pyclone_input.tsv
      python3 $moduleDir/pyclone_ctree.py --joint $outDir/pyclone_input.tsv --best_fit $best_fit --ctree_input $path_ctree

      """
}
