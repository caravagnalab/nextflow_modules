process PYCLONEVI {
    publishDir (
      params.publish_dir,
      mode: "copy"
    )
               
    input:

      tuple val(datasetID), val(patientID), val(sampleID), path(joint_table) // from the formatter output

    output:
      // tuple val(patientID), val(sampleID), path(path_ctree), emit: ctree_input
      tuple val(datasetID), val(patientID), val(sampleID), path(all_fits), emit: pyclone_all_fits
      tuple val(datasetID), val(patientID), val(sampleID), path(best_fit), emit: pyclone_best_fit
      // tuple val(patientID), val(sampleID), path(pyclone_joint), emit: pyclone_anno_joint
    script:
      def args = task.ext.args ?: ''
      def n_cluster_arg                    = args.n_cluster                     ?  "$args.n_cluster" : ""
      def density_arg                    = args.density                     ?  "$args.density" : ""
      def n_grid_point_arg                    = args.n_grid_point                     ?  "$args.n_grid_point" : ""
      def n_restarts_arg                    = args.n_restarts                     ?  "$args.n_restarts" : ""
      def mode                    = args.mode                     ?  "$args.mode" : ""

      if (mode == "singlesample") {
        outDir = "subclonal_deconvolution/pyclonevi/$datasetID/$patientID/$sampleID/"
        outDir_ctree = "$datasetID/$patientID/$sampleID/ctree"
        all_fits = "subclonal_deconvolution/pyclonevi/$datasetID/$patientID/$sampleID/all_fits.h5"
        best_fit = "subclonal_deconvolution/pyclonevi/$datasetID/$patientID/$sampleID/best_fit.txt"
        path_ctree = "$datasetID/$patientID/$sampleID/ctree/ctree_input_pyclonevi.csv"
        pyclone_joint = "subclonal_deconvolution/pyclonevi/$datasetID/$patientID/$sampleID/pyclone_joint.tsv"
      } else if (mode == "multisample"){
        sampleID = sampleID.join(' ')
        outDir = "subclonal_deconvolution/pyclonevi/$datasetID/$patientID"
        outDir_ctree = "$datasetID/$patientID/ctree"
        all_fits = "subclonal_deconvolution/pyclonevi/$datasetID/$patientID/all_fits.h5"
        best_fit = "subclonal_deconvolution/pyclonevi/$datasetID/$patientID/best_fit.txt"
        path_ctree = "$datasetID/$patientID/ctree/ctree_input_pyclonevi.csv"
        pyclone_joint = "subclonal_deconvolution/pyclonevi/$datasetID/$patientID/pyclone_joint.tsv"
      }

      """

      mkdir -p $outDir
      mkdir -p $outDir_ctree
      
      # format the input table in order to be pyclone compliant
      python3 $moduleDir/pyclone_utils.py create_pyclone_input $joint_table $patientID pyclone_input.tsv
    
      colnames="mutation_id\tsample_id\tref_counts\talt_counts\tnormal_cn\tmajor_cn\tminor_cn\ttumour_content"
      echo -e "\$colnames" > $outDir/pyclone_input.tsv 
      for i in $sampleID;
        do awk '\$2 == "'"\$i"'"' pyclone_input.tsv >> $outDir/pyclone_input.tsv;
      done
      
      pyclone-vi fit -i $outDir/pyclone_input.tsv -o $all_fits -c $n_cluster_arg -d $density_arg --num-grid-points $n_grid_point_arg --num-restarts $n_restarts_arg
      pyclone-vi write-results-file -i $all_fits -o $best_fit

      # python3 $moduleDir/pyclone_ctree.py --joint $outDir/joint_table.tsv --best_fit $best_fit --ctree_input $path_ctree --pyclone_joint $pyclone_joint

      """
}
