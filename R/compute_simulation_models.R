#' Train a model to simulate neurites.
#'
#' Train a model to simulate neurites. Learn two models, one to place the next node when the neurite continues and another one when the neurite bifurcates
#'
#' @param data is a dataframe which contains the features of each node
#' @param num_restart is an integer determining the number of restarts to learn the BN
#' @param max_iter maximum number of iterations to learn the BN
#'
#' @return a list with two models, the first one is the continuation model and the second one the bifurcation model
#'
#' @examples
#' path<-"/home/universidad/Documents/neuron_data/datos/All"
#' data<-get_features_files(path,60)
#' data<-compute_clusters(data)
#' simulation_model<-compute_simulation_models(data)
compute_simulation_models<-function(data,num_restart=5,max_iter=Inf)
{
  #Copy data and cast variables to factor
  cp_data<-data[,which(colnames(data)!="node_num_descendant")]
  cp_data$node_order<-as.factor(as.numeric(as.character(cp_data$node_order)))
  cp_data$desc_longer<-as.factor(as.numeric(as.character(cp_data$desc_longer)))
  cp_data$dendrite_diameter<-as.factor(as.numeric(as.character(cp_data$dendrite_diameter)))

  #Generate two data.frames, one for continue data and another for bifurcation data
  continue_data<-cp_data[which(data$node_num_descendant==1),]
  bifurcation_data<-cp_data[which(data$node_num_descendant==2),]

  #Arrays with the name of the variables to simulate
  first_desc<-c("desc_azimuth_angle","desc_elevation_angle","desc_length")
  second_desc<-c("desc_azimuth_angle2","desc_elevation_angle2","desc_length2")

  data_names<-colnames(cp_data)
  non_desc_var<-setdiff(data_names,c(first_desc,second_desc,"node_num_descendant","desc_longer","dendrite_diameter"))

  #Forbid arcs from continuous variables to factor variables
  factor_var<-intersect(data_names,c("node_order","subtree_max_order","subtree_min_order","subtree_no_bif","subtree_terminals","desc_longer","dendrite_diameter"))
  continuous_var<-setdiff(data_names,factor_var)
  black_list<-expand.grid(continuous_var,factor_var)
  black_list<-unique(rbind(black_list,expand.grid(first_desc,non_desc_var)))

  #Generate the BN model for continuation data
  continue_data<-continue_data[,!colnames(continue_data)%in%c("desc_azimuth_angle2","desc_elevation_angle2","desc_length2","desc_longer","dendrite_diameter")]
  black_list_c<-black_list[which(!black_list[,1]%in%c(second_desc,"desc_longer") & !black_list[,2]%in%c(second_desc,"desc_longer","dendrite_diameter")),]
  model_continue<-list()
  model_continue$structure<-bnlearn:::hc(continue_data,blacklist=black_list_c,restart=num_restart,max.iter = max_iter)
  model_continue$params<-bn.fit(model_continue$structure,continue_data)

  #Generate the BN model for bifurcation data
  black_list_b<-as.matrix(expand.grid(c(first_desc,second_desc),non_desc_var))
  black_list_b<-unique(rbind(black_list_b,expand.grid(continuous_var,factor_var)))
  model_bifurcation<-list()
  model_bifurcation$structure<-bnlearn:::hc(bifurcation_data,blacklist=black_list_b,restart=num_restart,max.iter = max_iter)
  model_bifurcation$params<-bn.fit(model_bifurcation$structure,bifurcation_data)

  return(list(continue=model_continue,bifurcation=model_bifurcation))
}

#' Dendrite diamater was estimate by this function. It is consider a latent variable in the model
#'
#' Dendrite diamater was estimate by this function. It is consider a latent variable in the model
#'
#' @param data is a data.frame with the measures computed in the terminal tips
#'
#' @return a data.frame that is the data input param after adding the new diameter column
compute_clusters<-function(data)
{
  temp_d <- data[(data$desc_length > 0) & (data$desc_length2 > 0),]
  cluster <- Mclust(temp_d$desc_length + temp_d$desc_length2)
  data$dendrite_diameter <- rep(1, nrow(data))
  data$dendrite_diameter[(data$desc_length > 0) & (data$desc_length2 > 0)] <- cluster$classification
  return(data)
}
