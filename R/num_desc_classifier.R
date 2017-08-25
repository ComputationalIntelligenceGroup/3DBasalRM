#' This function compute the features of each node
#'
#' This function compute the features of each node. Thus, for each node defining a neurite a set of features is computed according to node_feature_extractor function
#' @param path is a string indicating the folder where the files representing neurons are saved
#' @param eps is an integer value that represents the level of simplification applied over the neuron
#'
#' @return a dataset with the features of each node. List if features can be found in c_node_feature_extractor
#' @examples
#' path<-"/home/universidad/Documents/neuron_data/datos/All"
#' get_features_files(path,60)
get_features_files <- function(path, eps=60)
{
  factor_var <- c("node_num_descendant", "node_order", "subtree_max_order", "subtree_min_order", "subtree_no_bif", "subtree_terminals")

  path2files <- file.path(path, list.files(path))
  data <- data.table()

  for(file_path in path2files)
  {
    neuron <- neuro_converter(file_path, eps=eps)
    features <- rbindlist(node_feature_extractor(neuron))
    data <- rbind(data, features)
  }
  setDF(data)

  data$node_order[as.numeric(as.character(data$node_order)) > 4] <- 4
  for(i in factor_var)
  {
    data[, i] <- as.factor(data[, i])
  }

  #Remove segments of 0 length
  idx_zero <- which((as.numeric(as.character(data$node_num_descendant)) > 0) & (data$desc_length <= 0))
  if(length(idx_zero) > 0)
  {
    data <- data[-idx_zero, ]
  }
  data$desc_length <- log(data$desc_length)

  idx_zero <- which((as.numeric(as.character(data$node_num_descendant)) == 2) & (data$desc_length2 <= 0))
  if(length(idx_zero) > 0)
  {
    data<-data[-idx_zero, ]
  }

  data$desc_length2 <- log(data$desc_length2)

  return(data)
}


#' This function train a model to distinguish between those nodes that should be terminal, bifurcation or continue the neurite (0, 1 or 2 descendants)
#'
#' This function train a model to distinguish between those nodes that should be terminal, bifurcation or continue the neurite (0, 1 or 2 descendants)
#'
#' @param data is a dataframe which contains the features of each node
#'
#' @return model is an augmented Bayesian network classifier
#'
#' @examples
#' path<-"/home/universidad/Documents/neuron_data/datos/All"
#' data<-get_features_files(path,60)
#' desc_model<-train_num_descendant_model(data)
train_num_descendant_model <- function(data)
{
  data <- data[as.numeric(as.character(data$node_num_descendant))<3,]#Remove trifurcations
  attributes(data$node_num_descendant)$levels <- c("0", "1", "2")
  data <- data[as.numeric(as.character(data$node_order)) > 0, ]#Remove nodes with order 0 because they are deterministic, always continue.
  # data$node_order[as.numeric(as.character(data$node_order))>2]<-2 #Nodes with order bigger than 2 are grouped with 2
  featured_data<-data[,colnames(data) %in% c("azimuth_angle", "compartment_length", "distance_to_root", "elevation_angle", "length_to_brach_root", "node_num_descendant", "node_order", "node_to_brach_root_dist", "path_to_root", "subtree_box_volume", "subtree_length", "tortuosity")]#"path_to_root,subtree_box_volume"

  colnames(featured_data)[which(colnames(featured_data)=="node_num_descendant")] <- "class"
  featured_data$class <- factor(featured_data$class, levels=as.numeric(attributes(table(featured_data$class))$dimnames[[1]]))
  featured_data$node_order <- factor(featured_data$node_order, levels=as.numeric(attributes(table(featured_data$node_order))$dimnames[[1]]))

  model<-learn_BN(featured_data)

  return(model)
}
