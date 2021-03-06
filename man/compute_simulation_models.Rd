% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_simulation_models.R
\name{compute_simulation_models}
\alias{compute_simulation_models}
\title{Train a model to simulate neurites.}
\usage{
compute_simulation_models(data, num_restart = 5, max_iter = Inf)
}
\arguments{
\item{data}{is a dataframe which contains the features of each node}

\item{num_restart}{is an integer determining the number of restarts to learn the BN}

\item{max_iter}{maximum number of iterations to learn the BN}
}
\value{
a list with two models, the first one is the continuation model and the second one the bifurcation model
}
\description{
Train a model to simulate neurites. Learn two models, one to place the next node when the neurite continues and another one when the neurite bifurcates
}
\examples{
path<-"/home/universidad/Documents/neuron_data/datos/All"
data<-get_features_files(path,60)
data<-compute_clusters(data)
simulation_model<-compute_simulation_models(data)
}
