#' Get reconstruction from a file
#'
#' Given a reconstruction file, the converter guess the format by the file extension and reads it.
#' Optionally, it tries to correct errors in the reconstruction and applies a simplification algorithm over the branches.
#' Finally reconstruction is get in R.
#'
#' @param filepath a path to the file from features must be computed
#' @param correct a boolean value. The converter calls the corect method on each neuron in the reconstruction
#' @param eps is a boolean value. Error tolerance for the Ramer-Douglas-Peucker simplification algorithm applied at branch level.
#'
#' @return a data.frame containing the values measured for each neurite in the reconstruction
neuro_converter<-function(file_path, correct=F, eps=0.0)
{
  neuron_json <- c_neuro_converter(file_path, correct, eps)

  neuron <- fromJSON(neuron_json)

  return(list(plain=neuron_json, data=neuron))
}
