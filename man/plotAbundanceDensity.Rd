% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/plotAbundanceDensity.R
\name{plotAbundanceDensity}
\alias{plotAbundanceDensity}
\alias{plotAbundanceDensity,SummarizedExperiment-method}
\title{Plot abundance density}
\usage{
plotAbundanceDensity(x, ...)

\S4method{plotAbundanceDensity}{SummarizedExperiment}(
  x,
  layout = c("jitter", "density", "point"),
  assay.type = assay_name,
  assay_name = "counts",
  n = min(nrow(x), 25L),
  colour.by = colour_by,
  colour_by = NULL,
  shape.by = shape_by,
  shape_by = NULL,
  size.by = size_by,
  size_by = NULL,
  decreasing = order_descending,
  order_descending = TRUE,
  ...
)
}
\arguments{
\item{x}{a
\code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
object.}

\item{...}{additional parameters for plotting.
\itemize{
\item \code{xlab} \code{Character scalar}. Selects the x-axis label.
(Default: \code{assay.type})

\item \code{ylab} \code{Character scalar}. Selects the y-axis label.
\code{ylab} is disabled when \code{layout = "density"}.
(Default: \code{"Taxa"})

\item \code{point.alpha} \code{Numeric scalar}. From range 0 to 1. Selects
the transparency of
colour in \code{jitter} and \code{point} plot. (Default: \code{0.6})

\item \code{point.shape} \code{Positive integer scalar}. Value selecting
the shape of point in
\code{jitter} and \code{point} plot. (Default: \code{21})

\item \code{point.size} \code{Positive integer scalar}. Selects the size of
point in
\code{jitter} and \code{point} plot. (Default: \code{2})

\item \code{add_legend} \code{Logical scalar}. Determines if legend is
added. (Default: \code{TRUE})

\item \code{flipped}: \code{Logical scalar}. Determines if the orientation
of plot is changed so that x-axis and y-axis are swapped.
(Default: \code{FALSE})

\item \code{add_x_text} \code{Logical scalar}. Determines if text that
represents values is included in x-axis. (Default: \code{TRUE})

\item \code{jitter.height} \code{Numeric scalar}. Controls jitter in a
jitter plot. (Default: \code{0.25})

\item \code{jitter.width} \code{Numeric scalar}. Controls jitter in a
jitter plot. (Default: \code{NULL})
}
See \code{\link{mia-plot-args}} for more details i.e. call
\code{help("mia-plot-args")}}

\item{layout}{\code{Character scalar}. Selects the layout of the plot.
There are three different options: \code{jitter}, \code{density}, and
\code{point} plot. (default: \code{layout = "jitter"})}

\item{assay.type}{\code{Character scalar} value defining which assay data to
use. (Default: \code{"relabundance"})}

\item{assay_name}{Deprecate. Use \code{assay.type} instead.}

\item{n}{\code{Integer scalar}. Specifies the number of the most abundant
taxa to show. (Default: \code{min(nrow(x), 25L)})}

\item{colour.by}{\code{Character scalar}. Defines a column from
\code{colData}, that is used to color plot. Must be a value of
\code{colData()} function. (Default: \code{NULL})}

\item{colour_by}{Deprecated. Use \code{colour.by} instead.}

\item{shape.by}{\code{Character scalar}. Defines a column from
\code{colData}, that is used to group observations to different point shape
groups. Must be a value of \code{colData()} function. \code{shape.by} is
disabled when \code{layout = "density"}. (Default: \code{NULL})}

\item{shape_by}{Deprecated. Use \code{shape.by} instead.}

\item{size.by}{\code{Character scalar}. Defines a column from
\code{colData}, that is used to group observations to different point size
groups. Must be a value of \code{colData()} function. \code{size.by} is
disabled when \code{layout = "density"}. (Default: \code{NULL})}

\item{size_by}{Deprecated. Use \code{size.by} instead.}

\item{decreasing}{\code{Logical scalar}. Indicates whether the results should
be ordered in a descending order or not. If \code{NA} is given the order
as found in \code{x} for the \code{n} most abundant taxa is used.
(Default: \code{TRUE})}

\item{order_descending}{Deprecated. Use \code{order.descending} instead.}
}
\value{
A \code{ggplot2} object
}
\description{
This function plots abundance of the most abundant taxa.
}
\details{
This function plots abundance of the most abundant taxa. Abundance can be
plotted as a jitter plot, a density plot, or a point plot. By default, x-axis
represents abundance and y-axis taxa. In a jitter and point plot, each point
represents abundance of individual taxa in individual sample. Most common
abundances are shown as a higher density.

A density plot can be seen as a smoothened bar plot. It visualized
distribution of abundances where peaks represent most common abundances.
}
\examples{
data("peerj13075", package = "mia")
tse <- peerj13075

# Plots the abundances of 25 most abundant taxa. Jitter plot is the default
# option.
plotAbundanceDensity(tse, assay.type = "counts")

# Counts relative abundances
tse <- transformAssay(tse, method = "relabundance")

# Plots the relative abundance of 10 most abundant taxa.
# "nationality" information is used to color the points. X-axis is
# log-scaled.
plotAbundanceDensity(
    tse, layout = "jitter", assay.type = "relabundance", n = 10,
    colour.by = "Geographical_location") +
    scale_x_log10()

# Plots the relative abundance of 10 most abundant taxa as a density plot.
# X-axis is log-scaled
plotAbundanceDensity(
    tse, layout = "density", assay.type = "relabundance", n = 10 ) +
    scale_x_log10()

# Plots the relative abundance of 10 most abundant taxa as a point plot.
# Point shape is changed from default (21) to 41.
plotAbundanceDensity(
    tse, layout = "point", assay.type = "relabundance", n = 10,
    point.shape = 41)

# Plots the relative abundance of 10 most abundant taxa as a point plot.
# In addition to colour, groups can be visualized by size and shape in point
# plots, and adjusted for point size
plotAbundanceDensity(
    tse, layout = "point", assay.type = "relabundance", n = 10,
    shape.by = "Geographical_location", size.by = "Age", point.size=1)

# Ordering via decreasing
plotAbundanceDensity(
    tse, assay.type = "relabundance", decreasing = FALSE)

# for custom ordering set decreasing = NA and order the input object
# to your wishes
plotAbundanceDensity(
    tse, assay.type = "relabundance", decreasing = NA)

# Box plots and violin plots are supported by scater::plotExpression.
# Plots the relative abundance of 5 most abundant taxa as a violin plot.
library(scater)
top <- getTop(tse, top = 5)
plotExpression(tse, features = top, assay.type = "relabundance") +
    ggplot2::coord_flip()

# Plots the relative abundance of 5 most abundant taxa as a box plot.
plotExpression(tse, features = top, assay.type = "relabundance",
    show_violin = FALSE, show_box = TRUE) + ggplot2::coord_flip()

}
\seealso{
\code{\link[scater:plotExpression]{scater::plotExpression}}
}
\author{
Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
}
