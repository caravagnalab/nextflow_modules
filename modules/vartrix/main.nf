process VARTRIX {
    publishDir params.publish_dir, mode: 'copy'

    input:

      tuple val(patientID), val(sampleID), path(bam_File), path(bai_File), path(vcf_File), path(barcode_File)

    output:

      tuple path("$patientID/$sampleID/vartrix/*.mtx"), path("$patientID/$sampleID/vartrix/*.txt")

    script:
      def args = task.ext.args ?: ''
      def scoring_arg                    = args.scoring                     ?  "$args.scoring" : ""
      def padding_arg                    = args.padding                     ?  "$args.padding" : ""
      def umi_arg                    = args.umi                     ?  "$args.umi" : ""
      def refmatrix_arg = ""
      def output_prefix = "consensus"

      if (scoring_arg=='coverage'){
        refmatrix_arg = "${patientID}/${sampleID}/vartrix/ref_matrix_coverage.mtx"
        output_prefix = "coverage" 
      }

      """
      mkdir -p $patientID/$sampleID/vartrix


      vartrix -v $vcf_File \\
      -b $bam_File \\
      -f $params.ref_genome_vartrix \\
      -c $barcode_File \\
      -o $patientID/$sampleID/vartrix/matrix_${output_prefix}.mtx \\
      --out-variants $patientID/$sampleID/vartrix/variants_${output_prefix}.txt \\
      --scoring-method $scoring_arg \\
      --ref-matrix $refmatrix_arg \\
      --padding $padding_arg \\
      --umi $umi_arg
      """
}
