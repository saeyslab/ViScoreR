% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ViScoreR.R
\name{ComputeKNN}
\alias{ComputeKNN}
\title{Compute k-nearest-neighbour graph}
\usage{
ComputeKNN(X, k = 1000)
}
\arguments{
\item{X}{numeric matrix: row-wise coordinate matrix}

\item{k}{integer: size of neighbourhood (default: 1000)}
}
\description{
Finds k nearest neighbours to each point in a set of points using an approximate algorithm.
}
\references{
Lun A (2023). BiocNeighbors: Nearest Neighbor Detection for Bioconductor Packages. R package version 1.18.0.
}
