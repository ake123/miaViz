% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/plotCCA.R
\name{plotCCA}
\alias{plotCCA}
\alias{plotRDA}
\alias{plotCCA,SingleCellExperiment-method}
\alias{plotCCA,matrix-method}
\alias{plotRDA,SingleCellExperiment-method}
\alias{plotRDA,matrix-method}
\title{Plot RDA or CCA object}
\usage{
plotCCA(x, ...)

plotRDA(x, ...)

\S4method{plotCCA}{SingleCellExperiment}(x, dimred, ...)

\S4method{plotCCA}{matrix}(x, ...)

\S4method{plotRDA}{SingleCellExperiment}(x, dimred, ...)

\S4method{plotRDA}{matrix}(x, ...)
}
\arguments{
\item{x}{a
\code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-constructor]{TreeSummarizedExperiment}}
or a matrix of weights. The latter is returned as output from
\code{\link[mia:runCCA]{getRDA}}.}

\item{...}{additional parameters for plotting, inherited from
\code{\link[scater:plotReducedDim]{plotReducedDim}},
\code{\link[ggplot2:geom_label]{geom_label}} and
\code{\link[ggrepel:geom_label_repel]{geom_label_repel}}.
\itemize{
\item \code{add.ellipse}: One of
\code{c(TRUE, FALSE, "fill", "colour")}, indicating whether
ellipses should be present, absent, filled or colored.
(default: \code{ellipse.fill = TRUE})

\item \code{ellipse.alpha}: \code{Numeric scalar}. Between 0 and 1.
Adjusts the opacity of ellipses. (Default: \code{0.2})

\item \code{ellipse.linewidth}: \code{Numeric scalar}. Specifies the size
of ellipses. (Default: \code{0.1})

\item \code{ellipse.linetype}: \code{Integer scalar}. Specifies the style
of ellipses. (Default: \code{1})

\item \code{confidence.level}: \code{Numeric scalar}. Between 0 and 1.
Adjusts confidence level. (Default: \code{0.95})

\item \code{add.vectors}: \code{Logical scalar} or \code{character vector}.
If boolean, should vectors appear in the plot. If character,
selects vectors that are showed. The matching is done with regular
expression. (Default: \code{TRUE})

\item \code{vec.size}: \code{Numeric scalar}. Specifies the size of
vectors. (Default: \code{0.5})

\item \code{vec.colour}: \code{Character scalar}. Specifies the colour of
vectors. (Default: \code{"black"})

\item \code{vec.linetype}: \code{Integer scalar}. Specifies the style of
vector lines. (Default: \code{1})

\item \code{arrow.size}: \code{Numeric scalar}. Specifies the size of
arrows. (Default: \code{arrow.size = 0.25})

\item \code{label.size}: \code{Numeric scalar}. Specifies the size of text
and labels. (Default: \code{4})

\item \code{label.colour}: \code{Character scalar}. Specifies the colour of
text and labels. (Default: \code{"black"})

\item \code{sep.group}: \code{Character scalar}. Specifies the separator
used in the labels. (Default: \code{"\U2014"})

\item \code{repl.underscore}: \code{Character scalar}. Used to replace
underscores in the labels. (Default: \code{" "})

\item \code{vec.text}: \code{Logical scalar}. Should text instead of labels
be used to label vectors. (Default: \code{TRUE})

\item \code{repel.labels}: \code{Logical scalar}. Should labels be
repelled. (Default: \code{TRUE})

\item \code{parse.labels}: \code{Logical scalar}. Should labels be parsed.
(Default: \code{TRUE})

\item \code{add.significance}: \code{Logical scalar}. Should explained
variance and p-value appear in the labels. (Default: \code{TRUE})

\item \code{add.expl.var}: \code{Logical scalar}. Should explained
variance appear on the coordinate axes. (Default: \code{FALSE})

\item \code{add.centroids}: \code{Logical scalar}. Should centroids
of variables be added. (Default: \code{FALSE})

\item \code{add.species}: \code{Logical scalar}. Should species
scores be added. (Default: \code{FALSE})
}}

\item{dimred}{\code{Character scalar} or \code{integer scalar}. Determines
the reduced dimension to
plot. This is the output of \code{\link[mia:runCCA]{addRDA}} and resides in
\code{reducedDim(tse, dimred)}.}
}
\value{
A \code{ggplot2} object
}
\description{
\code{plotRDA} and \code{plotCCA} create an RDA/CCA plot starting from the
output of \code{\link[mia:runCCA]{CCA and RDA}} functions, two common methods
for supervised ordination of microbiome data.
}
\details{
\code{plotRDA} and \code{plotCCA} create an RDA/CCA plot starting from the
output of \code{\link[mia:runCCA]{CCA and RDA}} functions, two common methods
for supervised ordination of microbiome data. Either a
\code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-constructor]{TreeSummarizedExperiment}}
or a matrix object is supported as input. When the input is a
TreeSummarizedExperiment, this should contain the output of addRDA
in the reducedDim slot and the argument \code{dimred} needs to be defined.
When the input is a matrix, this should be returned as output from
getRDA. However, the first method is recommended because it provides
the option to adjust aesthetics to the colData variables through the
arguments inherited from \code{\link[scater:plotReducedDim]{plotReducedDim}}.
}
\examples{
# Load dataset
library(miaViz)
data("enterotype", package = "mia")
tse <- enterotype

# Run RDA and store results into TreeSE
tse <- transformAssay(tse, method = "relabundance")
tse <- addRDA(
    tse,
    assay.type = "relabundance",
    formula = assay ~ ClinicalStatus + Gender + Age,
    distance = "bray",
    na.action = na.exclude
    )

suppressWarnings({
# Create RDA plot coloured by variable
plotRDA(tse, "RDA", colour.by = "ClinicalStatus")

# Create RDA plot with empty ellipses
plotRDA(tse, "RDA", colour.by = "ClinicalStatus", add.ellipse = "colour")

# Create RDA plot with text encased in labels
plotRDA(tse, "RDA", colour.by = "ClinicalStatus", vec.text = FALSE)

# Create RDA plot without repelling text
plotRDA(tse, "RDA", colour.by = "ClinicalStatus", repel.labels = FALSE)

# Create RDA plot without vectors
plotRDA(tse, "RDA", colour.by = "ClinicalStatus", add.vectors = FALSE)

# Calculate RDA as a separate object
rda_mat <- getRDA(
    tse,
    assay.type = "relabundance",
    formula = assay ~ ClinicalStatus + Gender + Age,
    distance = "bray",
    na.action = na.exclude
    )

# Create RDA plot from RDA matrix
plotRDA(rda_mat)
})

}
