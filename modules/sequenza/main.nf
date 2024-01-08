process SEQUENZA_EXTRACT {
  publishDir params.publish_dir, mode: 'copy'

  input:
    
    tuple val(datasetID), val(patientID), val(sampleID), path(T_bamFile), path(T_baiFile), path(N_bamFile), path(N_baiFile)
  
  output:
  
    path("$datasetID/$patientID/$sampleID/SEQUENZA_EXTRACT/*.gz", "$datasetID/$patientID/$sampleID/SEQUENZA_EXTRACT/*.seqz")

  script:

  """  

  mkdir -p $params.publish_dir/$datasetID/$patientID/$sampleID/SEQUENZA_EXTRACT/

  sequenza-utils gc_wiggle -w 50 --fasta $params.ref_genome  -o Homo_sapiens_assembly38.gc50Base.wig

  sequenza-utils bam2seqz -n $N_bamFile -t $T_bamFile -F $params.ref_genome -gc Homo_sapiens_assembly38.gc50Base.wig -o tumor_vs_normal.seqz.gz

  sequenza-utils seqz_binning --seqz tumor_vs_normal.seqz.gz -w 50 -o tumor_vs_normal.binned.seqz

  """

}
