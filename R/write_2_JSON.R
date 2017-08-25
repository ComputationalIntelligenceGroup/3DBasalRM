#' Write the neuron into a JSON file
#'
#' This function gets a neuron an writes it into a JSON file
#'
#' @param neuron info about a reconstruction in JSON format
#' @param file_path is a string with the path where the json file must be saved, it must include the name of the file and the extension
#'
#' @return NULL
#'
#' @export
write_2_JSON <- function(neuron, file_path="~/datos/simpl.json")
{
  write(neuron$plain, file_path)
}
