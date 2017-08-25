#' Main function. Execute the repairing process
#'
#' Main function. Execute the repairing process
#'
#' @param file_path is a string representing the path to the .DAT or JSON file containing Neurolucida data
#' @param seed is an integer value denoting the seed used to generate random variables
#'
#' @examples
#' file_path <- "/home/universidad/datos/Cortadas/human cing id if6 porta 2 sec1 cel20.DAT"
#' file_path <- "/home/universidad/datos/H213III/h213III1.DAT"
#' repaired_neuron <- main(file_path, 1)
#' write_2_JSON(repaired_neuron,path="~/datos/repaired_neuron.json")
#'
#' @export
main<-function(file_path, seed = 1)
{
  max_id_node <- 0; #Maximo id de los nodos que pertenecen a la neurona. Necesario para generar los nuevos ids de los nodos.

  set.seed(seed)

  #Load models
  desc_model <- neurostr::desc_model
  simulation_model <- neurostr::simulation_model

  #Read neuron
  neuron <- neuro_converter(file_path, eps=60)

  #Identify cut nodes and get their ids
  id_cut_nodes <- get_cut_nodes(neuron, 5)

  #Repeat until the list of cut nodes is empty
  while(length(id_cut_nodes) > 0)
  {
    features <- rbindlist(node_feature_extractor(neuron, id_nodes=id_cut_nodes))

    selected_features <- features[ ,colnames(features) %in% names(desc_model$params), with=F]
    selected_features$node_order[selected_features$node_order > 4] <- 4
    selected_features$node_order <- factor(selected_features$node_order, levels=0:4)


    #Classify nodes among bifurcation, continuation or terminal
    desc_probabilities <- pred_BN(desc_model, selected_features)

    #Sampled number of descendant radonmly according to the probability of terminal, continue o bifurcation for each node
    num_of_descendant <- apply(desc_probabilities[, grep("prob", colnames(desc_probabilities))], 1, function(x) {sample(c(0, 1, 2), size=1, prob=x)})

    #If there is at least one dendrite that is not terminal, make it grow
    if(sum(num_of_descendant) > 0)
    {
      #Simulate new nodes
      simulation_data <- features[, colnames(features) %in% names(simulation_model$continue$structure$nodes), with=F]
      simulation_data$node_order[simulation_data$node_order > 4] <- 4
      simulation_data$node_order <- factor(simulation_data$node_order, levels=0:4)

      #If there is at least one continuation save the information of the new node to add at the end of the dendrite
      if(nrow(simulation_data[num_of_descendant==1]) > 0)
      {
        continue_data <- simulation_data[num_of_descendant==1]
        simulation_continue <- simulate_continuation(simulation_model$continue,continue_data)
        simulation_continue$id_parent <- id_cut_nodes[num_of_descendant==1]
        simulated_data <- simulation_continue
      }
      #If there is at least one bifurcation save the information of the position of each new dendrite
      if(nrow(simulation_data[num_of_descendant==2]) > 0)
      {
        bifurcation_data <- simulation_data[num_of_descendant==2]
        simulation_bifurcation <- simulate_bifurcation(simulation_model$bifurcation,bifurcation_data)
        second_branch_idx <- grep("2", colnames(simulation_bifurcation))
        simulation_bifurcation_2 <- simulation_bifurcation[,second_branch_idx]
        colnames(simulation_bifurcation_2) <- gsub("2", "", colnames(simulation_bifurcation_2))
        simulation_bifurcation <- rbind(simulation_bifurcation[, setdiff(1:ncol(simulation_bifurcation), second_branch_idx)], simulation_bifurcation_2)
        simulation_bifurcation$id_parent <- rep(id_cut_nodes[num_of_descendant==2], 2)

        if(nrow(simulation_data[num_of_descendant==1]) > 0)
        {
          simulated_data <- rbind(simulation_continue, simulation_bifurcation)
        }else{
          simulated_data <- simulation_bifurcation
        }

      }

      #Recover original values of the length (it has been modeled as a log-norm)
      simulated_data$desc_length <- exp(simulated_data$desc_length)
      is_bifurcation <- c(rep(1, nrow(simulation_data[num_of_descendant==1])), rep(2, nrow(simulation_data[num_of_descendant==2]) * 2)) - 1


      #Add new nodes to the cut dendrites
      neuron <- simulated_node_coordinates(neuron$plain, data.matrix(simulated_data), is_bifurcation, max_id_node)
      id_cut_nodes <- neuron$ids

      max_id_node <- max(id_cut_nodes)
    }else{#If there was not a node to simulate the list of nodes is empty
      id_cut_nodes <- c()
    }
  }#End while
  return(neuron)
}
