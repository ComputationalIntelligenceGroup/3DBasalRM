#' Given a terminal node, simulate the grow of the dendrite one step
#'
#' Given a terminal node, simulate the grow of the dendrite one step. A step is a bifurcation or a continuation.
#'
#' @param model is a bn model generate by bnlearn
#' @param node_info is a data.frame containing the terminal nodes that must grow
#'
#' @return a data.frame with the spherical coordinates of the new node with respect to the terminal node introduced to the function
#'
#' @examples
#' path <- "/home/universidad/Documents/neuron_data/datos/All"
#' data <- get_features_files(path, 60)
#' data <- compute_clusters(data)
#' simulation_model <- compute_simulation_models(data)
#' simulation <- simulate_continuation(simulation_model$continue, data[data$node_num_descendant==1])
#'
#' @export
simulate_continuation <- function(model, node_info)
{
  #Get dataset of nodes to grow and get topological order of the variables according to the BN structure
  node_info <- data.frame(node_info)
  node_order_model <- node.ordering(model$structure)
  first_desc <- c("desc_azimuth_angle", "desc_elevation_angle", "desc_length")
  first_desc <- node_order_model[node_order_model %in% first_desc]

  #Save position of simulated nodes
  simulated_nodes <- matrix(0, ncol=length(first_desc), nrow=nrow(node_info))

  #For each feature to simulate
  for(i in 1:length(first_desc))
  {
    #Initiate mean and sd matrices to sample
    mean_matrix <- matrix(0, ncol=1, nrow=nrow(node_info))
    sd_matrix <- matrix(0, ncol=1, nrow=nrow(node_info))

    #Get the information of the parents of the node (are discrete and their names)
    node <- first_desc[i]
    parent_names <- parents(model$structure, node)
    is_discrete_p <- sapply(parent_names, function(x) {class(model$params[[x]])=="bn.fit.dnode"})
    discrete_p_names <- parent_names[is_discrete_p]
    continuous_p_names <- parent_names[!is_discrete_p]

    #If the node is Gaussian, that is, all its parents are Gaussian or dont have any parent
    if(class(model$params[[node]])=="bn.fit.gnode"){
      mean_matrix[,1] <- coefficients(model$params[[node]]) %*% t(cbind(1,node_info[,colnames(node_info) %in% continuous_p_names]))
      sd_matrix[,1] <- model$params[[node]]$sd
    }

    #If the node have discrete parents
    if(class(model$params[[node]])=="bn.fit.cgnode"){
      #Get the configuration of the parents to select the correct coefficients for the normal distribution
      pos_parents <- lapply(discrete_p_names, function(x) {match(node_info[,x], as.numeric(model$params[[node]]$dlevels[[x]]))})
      num_poss_parents <- sapply(model$params[[node]]$dlevels,length)#Number of possible values for each parent

      #The position in dlevels is a+num_poss_parents[1]*(b-1)+num_poss_parents[2]*num_poss_parents[1]*(c-1)+...
      match_idx <- pos_parents[[1]]#If there is only 1 parent then the position is a
      if(length(num_poss_parents) > 1)#If there is more than one parent the position is a combination
      {
        prev_parent_num_combinations <- 1
        for(j in 2:length(num_poss_parents))
        {
          prev_parent_num_combinations <- prev_parent_num_combinations * num_poss_parents[j-1]
          match_idx <- match_idx + prev_parent_num_combinations * (pos_parents[[j]]-1)
        }
      }

      #If all the parent nodes are discrete we only have the mean of the node
      if(length(parent_names)==length(discrete_p_names))
      {
        mean_matrix <- model$params[[node]]$coefficients[match_idx]
      }else{#If there is at least one gaussian parent the mean is a linear combination
        mean_matrix <- colSums(model$params[[node]]$coefficients[,match_idx] * t(cbind(1, node_info[,colnames(node_info) %in% continuous_p_names])))
      }
        sd_matrix <- model$params[[node]]$sd[match_idx]
    }

    #Simulate values for the node
    if(node=="desc_length")
    {
      simulated_nodes[,i] <- rtnorm(length(mean_matrix), mean=mean_matrix, sd=sd_matrix, upper=5.75)
    }else if(node=="desc_azimuth_angle")
      {
        simulated_nodes[,i] <- rtnorm(length(mean_matrix), mean=mean_matrix, sd=sd_matrix, lower=-pi, upper=pi)
      }else{
        simulated_nodes[,i] <- rtnorm(length(mean_matrix), mean=mean_matrix, sd=sd_matrix, lower=0, upper=pi)
      }

    #Save the simulation to use them for the next simulations
    node_info[[first_desc[i]]] <- simulated_nodes[,i]
    node_info <- node_info[,order(names(node_info))]
  }

  simulated_nodes <- data.frame(simulated_nodes)
  colnames(simulated_nodes) <- first_desc

  return(simulated_nodes)
}

#' Given a terminal node, simulate the bifurcation of the dendrite
#'
#' Given a terminal node, simulate the bifurcation of the dendrite
#'
#' @param model is a bn model generate by bnlearn
#' @param node_info is a data.frame containing the terminal nodes that must grow
#'
#' @return a data.frame with the spherical coordinates of the two new nodes with respect to the terminal node introduced to the function
#'
#' @examples
#' path <- "/home/universidad/Documents/neuron_data/datos/All"
#' data <- get_features_files(path, 60)
#' data <- compute_clusters(data)
#' simulation_model <- compute_simulation_models(data)
#' simulation <- simulate_bifurcation(simulation_model$bifurcation, data[data$node_num_descendant==2,])
#'
#' @export
simulate_bifurcation <- function(model, node_info)
{
  node_info <- data.frame(node_info)
  node_order_model <- node.ordering(model$structure)
  first_desc <- c("desc_azimuth_angle", "desc_elevation_angle", "desc_length", "desc_azimuth_angle2", "desc_elevation_angle2", "desc_length2")
  first_desc <- node_order_model[node_order_model %in% first_desc]

  #Simulate the discrete nodes
  disc_sim <- c("desc_longer", "dendrite_diameter")
  disc_sim <- node_order_model[node_order_model %in% disc_sim]
  for(i in 1:length(disc_sim))
  {
  node <- disc_sim[i]
  parents_longer <- parents(model$structure, node)
  if(length(parents_longer)==0)
  {
    node_info[[node]] <- c(rMultinom(matrix(model$params[[node]]$prob, nrow=1), nrow(node_info)))
  }else{
    #Get the configuration of the parents to select the correct probability
    pos_parents <- lapply(parents_longer, function(x) {match(node_info[,x], as.numeric(attributes(model$params[[node]]$prob)$dimnames[[x]]))})
    pos_parents <- c(list(rep(1, length(pos_parents[[1]]))), pos_parents)
    num_poss_parents <- sapply(attributes(model$params[[node]]$prob)$dimnames, length)#Number of possible values for each parent

    #The position in dlevels is a+num_poss_parents[1]*(b-1)+num_poss_parents[2]*num_poss_parents[1]*(c-1)+...
    match_idx <- pos_parents[[1]]#If there is only 1 parent then the position is a
    if(length(num_poss_parents) > 1)#If there is more than one parent the position is a combination
    {
      prev_parent_num_combinations <- 1
      for(j in 2:length(num_poss_parents))
      {
        prev_parent_num_combinations <- prev_parent_num_combinations * num_poss_parents[j-1]
        match_idx <- match_idx + prev_parent_num_combinations * (pos_parents[[j]]-1)
      }
    }
    match_matrix <- sapply(match_idx, function(x) {x:(x + num_poss_parents[1] - 1)})
    node_info[[node]] <- c(rMultinom(probs=matrix(model$params[[node]]$prob[c(match_matrix)], ncol=num_poss_parents[1],byrow=T), m=1))
  }
  node_info <- node_info[,order(names(node_info))]
  }
  #Save position of simulated nodes
  simulated_nodes <- matrix(0, ncol=length(first_desc), nrow=nrow(node_info))
  node_info <- node_info[,order(names(node_info))]

  #For each feature to simulate
  for(i in 1:length(first_desc))
  {
    #Initiate mean and sd matrices to save
    mean_matrix <- matrix(0, ncol=1, nrow=nrow(node_info))
    sd_matrix <- matrix(0, ncol=1, nrow=nrow(node_info))

    #Get the information of the parents of the node (are discrete and their names)
    node <- first_desc[i]
    parent_names <- parents(model$structure, node)
    is_discrete_p <- sapply(parent_names, function(x) {class(model$params[[x]])=="bn.fit.dnode"})
    discrete_p_names <- parent_names[is_discrete_p]
    continuous_p_names <- parent_names[!is_discrete_p]

    #If the node is Gaussian, that is, all its parents are Gaussian or dont have any parent
    if(class(model$params[[node]])=="bn.fit.gnode"){
      mean_matrix[,1] <- coefficients(model$params[[node]]) %*% t(cbind(1, node_info[,colnames(node_info) %in% continuous_p_names]))
      sd_matrix[,1] <- model$params[[node]]$sd
    }

    #If the node have discrete parents
    if(class(model$params[[node]])=="bn.fit.cgnode"){
      #Get the configuration of the parents to select the correct coefficients for the normal distribution
      pos_parents <- lapply(discrete_p_names, function(x) {match(node_info[,x], as.numeric(model$params[[node]]$dlevels[[x]]))})
      num_poss_parents <- sapply(model$params[[node]]$dlevels, length)#Numbre of possible values for each parent

      #The position in dlevels is a+num_poss_parents[1]*(b-1)+num_poss_parents[2]*num_poss_parents[1]*(c-1)+...
      match_idx <- pos_parents[[1]]#If there is only 1 parent then the position is a
      if(length(num_poss_parents) > 1)#If there is more than one parent the position is a combination
      {
        prev_parent_num_combinations <- 1
        for(j in 2:length(num_poss_parents))
        {
          prev_parent_num_combinations <- prev_parent_num_combinations * num_poss_parents[j-1]
          match_idx <- match_idx+prev_parent_num_combinations * (pos_parents[[j]] - 1)
        }
      }

      if(length(parent_names)==length(discrete_p_names))#If all the parent nodes are discrete we only have the mean of the node
      {
        mean_matrix <- model$params[[node]]$coefficients[match_idx]
      }else{#If there is at least one gaussian parent the mean is a linear combination
        mean_matrix <- colSums(model$params[[node]]$coefficients[,match_idx] * t(cbind(1,node_info[,colnames(node_info) %in% continuous_p_names])))
      }
      sd_matrix <- model$params[[node]]$sd[match_idx]
    }

    if(grepl("desc_length*", node)){
      simulated_nodes[,i] <- rtnorm(length(mean_matrix), mean=mean_matrix, sd=sd_matrix, upper=5.75)
    }else if(grepl("desc_azimuth_angle*", node)){
      simulated_nodes[,i] <- rtnorm(length(mean_matrix), mean=mean_matrix, sd=sd_matrix, lower=-pi, upper=pi)
    }else{
      simulated_nodes[,i] <- rtnorm(length(mean_matrix), mean=mean_matrix, sd=sd_matrix, lower=0, upper=pi)
    }
    node_info[[first_desc[i]]] <- simulated_nodes[,i]
    node_info <- node_info[,order(names(node_info))]
  }
  simulated_nodes <- data.frame(simulated_nodes)
  colnames(simulated_nodes) <- first_desc
  simulated_nodes <- simulated_nodes[,order(names(simulated_nodes))]

  return(simulated_nodes)
}

