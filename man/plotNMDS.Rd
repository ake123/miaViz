% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotNDMS.R
\name{plotNMDS}
\alias{plotNMDS}
\title{Wrapper for scater::plotReducedDim()}
\usage{
plotNMDS(x, ..., ncomponents = 2)
}
\arguments{
\item{x}{a
\code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
object.}

\item{...}{additional arguments passed to scater::plotReducedDim().}

\item{ncomponents}{\code{Numeric scalar}. indicating the number of dimensions
to plot, starting from the first dimension. Alternatively, a numeric vector
specifying the dimensions to be plotted. (Default: \code{2})}
}
\description{
Wrapper for scater::plotReducedDim()
}
