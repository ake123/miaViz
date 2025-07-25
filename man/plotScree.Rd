% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/plotScree.R
\name{plotScree}
\alias{plotScree}
\alias{plotScree,SingleCellExperiment-method}
\alias{plotScree,ANY-method}
\title{Create a scree plot}
\usage{
plotScree(x, ...)

\S4method{plotScree}{SingleCellExperiment}(x, dimred, ...)

\S4method{plotScree}{ANY}(x, ...)
}
\arguments{
\item{x}{a
\code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-constructor]{TreeSummarizedExperiment}}
\code{\link[vegan:eigenvals]{eigenvals}} or a vector.}

\item{...}{additional parameters for plotting
\itemize{
\item \code{show.barplot}: \code{Logical scalar}. Whether to show a
barplot. (Default: \code{TRUE})

\item \code{show.points}: \code{Logical scalar}. Whether to show a
points. (Default: \code{TRUE})

\item \code{show.line}: \code{Logical scalar}. Whether to show a
line. (Default: \code{TRUE})

\item \code{show.labels}: \code{Logical scalar}. Whether to show a
labels for each point. (Default: \code{FALSE})

\item \code{add.proportion}: \code{Logical scalar}. Whether to show
proportion of explained variance, i.e., raw eigenvalues.
(Default: \code{TRUE})

\item \code{add.cumulative}: \code{Logical scalar}. Whether to show
cumulative explained variance calculated from eigenvalues.
(Default: \code{FALSE})

\item \code{n}: \code{Integer scalar}. Number of eigenvalues to plot.
If \code{NULL}, all eigenvalues are plotted. (Default: \code{NULL})

\item \code{show.names}: \code{Logical scalar}. Whether to show names of
components in x-axis. If \code{FALSE}, the index of component is shown
instead of names. (Default: \code{FALSE})

\item \code{eig.name}: \code{Character scalar}. The name of the attribute
in \code{reducedDim(x, dimred)} that contains the eigenvalues.
(Default: \code{c("eig", "varExplained")})
}}

\item{dimred}{\code{Character scalar} or \code{integer scalar}. Determines
the reduced dimension to plot. This is used when \code{x} is a
\code{TreeSummarizedExperiment} to extract the eigenvalues from
\code{reducedDim(x, dimred)}.}
}
\value{
A \code{ggplot2} object
}
\description{
\code{plotScree} generates a scree plot to visualize the eigenvalues.
The eigenvalues can be provided either as a part of a
\code{TreeSummarizedExperiment} object or as a separate \code{vector}.
This plot illustrates the decline in eigenvalues across components,
helping to assess the importance of each component.
}
\details{
\code{plotScree} generates a scree plot to visualize the relative importance
of components in dimensionality reduction techniques such as Principal
Component Analysis (PCA) or Principal Coordinate Analysis (PCoA). If the
input is a \code{TreeSummarizedExperiment} object, the function extracts
eigenvalues from the specified reduced dimension slot, which requires that
dimensionality reduction has been performed beforehand using a dedicated
function. Alternatively, if the input is a \code{vector} or an
\code{eigenvals} object, these values are directly used as eigenvalues for
the plot.

The plot can include a combination of barplot, points, connecting lines,
and labels, which can be controlled using the \code{show.*} parameters.

An option to show cumulative explained variance is also available by setting
\code{add.cumulative = TRUE}.
}
\examples{

library(miaViz)
library(scater)

data("enterotype", package = "mia")
tse <- enterotype

# Run  PCA and store results into TreeSE
tse <- transformAssay(tse, method = "clr", pseudocount = TRUE)
tse <- runPCA(tse, assay.type = "clr")

# Plot scree plot
plotScree(tse, "PCA", add.cumulative = TRUE)

}
