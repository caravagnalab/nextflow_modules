process JOINT_FIT {
  publishDir (
    params.publish_dir,
    mode: "copy"
  )

  input:
    tuple val(datasetID), val(patientID), val(sampleID), path(mobster_best_fits), path(joint_table)

  output:
    tuple val(datasetID), val(patientID), val(sampleID), path("$outDir/mCNAqc_filtered.rds"), emit: annotated_joint

  script:
    def args = task.ext.args ?: ""

    outDir = "subclonal_deconvolution/mobster/$datasetID/$patientID"

    if (!(sampleID instanceof String)) {
      sampleID = sampleID.join(",")
    }

    if (!(mobster_best_fits instanceof String)) {
      mobster_best_fits = mobster_best_fits.join(",")
    }

    """
    #!/usr/bin/env Rscript

    # Sys.setenv("VROOM_CONNECTION_SIZE"=99999999)

    library(dplyr)

    patientID = "$patientID"
    # samples = strsplit(x="$sampleID", ",") %>% unlist()  # list of samples
    fits = strsplit(x="$mobster_best_fits", ",") %>% unlist()  # list of mobster fitnames

    muts_to_discard = lappply(fits, function(fit_name) {
      fit_i = readRDS(fit_name)
      fit_i[["data"]] %>% 
        dplyr::filter(cluster=="Tail") %>% 
        dplyr::mutate(sample_id=fit_name) %>% 
        dplyr::select(mutation_id, cluster, sample_id)
    }) %>% dplyr::bind_rows() %>% 
      dplyr::group_by(mutation_id) %>% 
      dplyr::summarise(n_samples=dplyr::n()) %>% 
      dplyr::filter(n_samples == length(fits)) %>%
      dplyr::pull(mutation_id)

    # Function to discard mutations from the mCNAqc objects
    discard_mutations_from_mCNAqc = function(obj, mutations_list) {
      sample_names = obj[["original_cnaqc_objc"]] %>% names()
      
      get_mCNAqc_types = function(x) {
        names(x) %>% purrr::discard(function(i) i == "m_cnaqc_stats")
      }
      
      types = get_mCNAqc_types(obj)
      
      for (tid in types) {
        sample_names = names(obj[[tid]])
        for (sample_id in sample_names) {
          obj[[tid]][[sample_id]][["mutations"]] = CNAqc::Mutations(obj[[tid]][[sample_id]]) %>% 
            # dplyr::mutate(unique_id=paste(chr, from, to, ref, alt, sep="_")) %>% 
            dplyr::filter(!mutation_id %in% mutations_list)
        }
      }
      
      return(obj)
    }

    if ( grepl(".rds\$", tolower("$joint_table")) ) {
      obj = readRDS("$joint_table")
      if (class(obj) == "m_cnaqc") {
        obj_filtered = discard_mutations_from_mCNAqc(obj, muts_to_discard)

        saveRDS(obj_filtered, "$outDir/mCNAqc_filtered.rds")
      } else {
        cli::cli_alert_warning("Object of class {class($joint_table)} not supported. Saving the original object.")
        saveRDS(obj, "$outDir/mCNAqc_filtered.rds")
      }
    }

    """
}
