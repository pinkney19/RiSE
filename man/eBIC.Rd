% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{eBIC}
\alias{eBIC}
\title{Function to evaluate eBIC}
\usage{
eBIC(theta, S, z, gam, n.stream, f, n.trials)
}
\arguments{
\item{theta}{estimate of inverse SDM for given lambda}

\item{S}{Peirodogram estimate}

\item{gam}{tuning parameter}

\item{n.stream}{number of dimensions/neurons}

\item{f}{number of frequencies considered in trial-frequency smoothed estimator}

\item{n.trials}{number of experimental trials (m)}
}
\description{
Function to evaluate eBIC
}
