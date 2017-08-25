#' Include the simulated nodes into the neuron
#'
#' Include the simulated nodes into the neuron representation.
#'
#' @param json_info info about a reconstruction in JSON format
#' @param simualted_nodes is a matrix with the spherical coordinates of the new nodes with respect to their parent nodes
#' @param is_bifurcation is a vector identifying if the node to simulate is a bifurcation or is a continuation
#' @param init_id is an integer denoting the highest id of the neuron, it is used to give ids to the new nodes
#' @param omit_apical a boolean value. If set, apical dendrite is not measured
#' @param omit_axon a boolean value. If set, axon is not measured
#' @param omit_dend a boolean value. If set, dendrites are not measured
#' @param correct a boolean value. The converter calls the correct method on each neuron in the reconstruction
#' @param removeZjumps a boolean value. Remove jumps over Z axis produced by the Z jumps in the microscopy imaging
#'
#' @return a list of three elements, the first one is the neuron in JSON format, the second one is a data.frame representing the neuron information and the third one are the ids given to the simulated nodes
#'
#' @export
simulated_node_coordinates <- function(json_info, simulated_nodes, is_bifurcation, init_id, omit_apical = F, omit_axon = F, omit_dend = F, correct = F, removeZjumps = T)
{
  neuron_json <- c_simulated_position(json_info, simulated_nodes, is_bifurcation, init_id, omit_apical, omit_axon, omit_dend, correct, removeZjumps)
  neuron <- fromJSON(neuron_json[[1]])

  return(list(plain = neuron_json[[1]], data = neuron, ids = neuron_json[[2]]))
}
