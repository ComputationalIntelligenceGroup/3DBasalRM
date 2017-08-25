#' Compute neurite features
#' The neurite feature extractor computes a set of prebuilt measure for each neurite in the reconstruction.
#' Measure values are not aggregated in any way, they are returned as-is (usually as vectors).
#' Optionally, it tries to correct errors in the reconstruction.
#'
#' @param json_info info about a reconstruction in JSON format
#' @param omit_apical a boolean value. If set, apical dendrite is not measured
#' @param omit_axon a boolean value. If set, axon is not measured
#' @param omit_dend a boolean value. If set, dendrites are not measured
#' @param correct a boolean value. The converter calls the corect method on each neuron in the reconstruction
#'
#' @return a data.frame containing the values measured for each neurite in the reconstruction
#'
#' @export
neurite_feature_extractor<-function(json_info, omit_apical=F, omit_axon=F, omit_dend=F, correct=F)
{
  neurite_features_raw <- c_neurite_feature_extractor(json_info, omit_apical, omit_axon, omit_dend, correct)

  neurite_features <- fromJSON(neurite_features_raw)

  return(neurite_features)
}
