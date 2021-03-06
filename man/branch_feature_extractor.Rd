% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/branch_feature_extractor.R
\name{branch_feature_extractor}
\alias{branch_feature_extractor}
\title{Compute neurite features
The neurite feature extractor computes a set of prebuilt measure for each brench in the reconstruction.
Optionally, it tries to correct errors in the reconstruction.}
\usage{
branch_feature_extractor(json_info, omit_apical = F, omit_axon = F,
  omit_dend = F, correct = F)
}
\arguments{
\item{json_info}{info about a reconstruction in JSON format}

\item{omit_apical}{a boolean value. If set, apical dendrite is not measured}

\item{omit_axon}{a boolean value. If set, axon is not measured}

\item{omit_dend}{a boolean value. If set, dendrites are not measured}

\item{correct}{a boolean value. The converter calls the corect method on each neuron in the reconstruction}
}
\value{
a data.frame containing the values measured for each neurite in the reconstruction
}
\description{
Compute neurite features
The neurite feature extractor computes a set of prebuilt measure for each brench in the reconstruction.
Optionally, it tries to correct errors in the reconstruction.
}
