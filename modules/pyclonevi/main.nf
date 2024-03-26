process PYCLONEVI {
    publishDir (
      params.publish_dir,
      mode: "copy"
    )
               
    input:

      tuple val(patientID), val(sampleID), path(joint_table)

    output:
      // tuple val(patientID), val(sampleID), path(path_ctree), emit: ctree_input
      tuple val(patientID), val(sampleID), path(all_fits), emit: pyclone_all_fits
      tuple val(patientID), val(sampleID), path(best_fit), emit: pyclone_best_fit
      tuple val(patientID), val(sampleID), path(pyclone_joint), emit: pyclone_anno_joint
    script:
      def args = task.ext.args ?: ''
      def n_cluster_arg                    = args.n_cluster                     ?  "$args.n_cluster" : ""
      def density_arg                    = args.density                     ?  "$args.density" : ""
      def n_grid_point_arg                    = args.n_grid_point                     ?  "$args.n_grid_point" : ""
      def n_restarts_arg                    = args.n_restarts                     ?  "$args.n_restarts" : ""
      def step                    = args.step                     ?  "$args.step" : ""

      if (step == "subclonal_singlesample") {
        outDir = "$patientID/$sampleID/pyclonevi"
        outDir_ctree = "$patientID/$sampleID/ctree"
        all_fits = "$patientID/$sampleID/pyclonevi/all_fits.h5"
        best_fit = "$patientID/$sampleID/pyclonevi/best_fit.txt"
        path_ctree = "$patientID/$sampleID/ctree/ctree_input_pyclonevi.csv"
        pyclone_joint = "$patientID/$sampleID/pyclonevi/pyclone_joint.tsv"
      } else if (step == "subclonal_multisample"){
        sampleID = sampleID.join(' ')
        outDir = "$patientID/pyclonevi"
        outDir_ctree = "$patientID/ctree"
        all_fits = "$patientID/pyclonevi/all_fits.h5"
        best_fit = "$patientID/pyclonevi/best_fit.txt"
        path_ctree = "$patientID/ctree/ctree_input_pyclonevi.csv"
        pyclone_joint = "$patientID/pyclonevi/pyclone_joint.tsv"
      }

      """

      mkdir -p $outDir
      mkdir -p $outDir_ctree
      echo $sampleID
      head -1 $joint_table | sed 's/NV/alt_counts/g' | sed 's/NR/ref_counts/g' | sed 's/purity/tumour_content/g' > $outDir/joint_table.tsv
      # head -1 $joint_table | sed 's/NV/alt_counts/g' | sed 's/NR/ref_counts/g' | sed 's/purity/tumour_content/g' | sed 's/patient_id/patientID/g' | sed 's/gene/variantID/g' | sed 's/is_driver/is.driver/g' | sed 's/driver_label/variantID/g'> $outDir/joint_table.tsv 
      for i in $sampleID;
      do awk '\$2 == "'"\$i"'"' $joint_table >> $outDir/joint_table.tsv;
      done

      pyclone-vi fit -i $outDir/joint_table.tsv -o $all_fits -c $n_cluster_arg -d $density_arg --num-grid-points $n_grid_point_arg --num-restarts $n_restarts_arg
      pyclone-vi write-results-file -i $all_fits -o $best_fit

      python3 $moduleDir/pyclone_ctree.py --joint $outDir/joint_table.tsv --best_fit $best_fit --ctree_input $path_ctree --pyclone_joint $pyclone_joint

      """
}
