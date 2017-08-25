#' Compute neurite features
#' The neurite feature extractor computes a set of prebuilt measure for each brench in the reconstruction.
#' Optionally, it tries to correct errors in the reconstruction.
#'
#' @param json_info info about a reconstruction in JSON format
#' @param omit_apical a boolean value. If set, apical dendrite is not measured
#' @param omit_axon a boolean value. If set, axon is not measured
#' @param omit_dend a boolean value. If set, dendrites are not measured
#' @param correct a boolean value. The converter calls the corect method on each neuron in the reconstruction
#'
#' @return a data.frame containing the values measured for each neurite in the reconstruction
branch_feature_extractor<-function(json_info, omit_apical=F, omit_axon=F, omit_dend=F, correct=F)
{
  branch_features_raw <- c_branch_feature_extractor(json_info, omit_apical, omit_axon, omit_dend, correct)

  branch_features <- fromJSON(branch_features_raw)

  return(branch_features)
}
