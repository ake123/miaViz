% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotGraph.R
\name{plotGraph}
\alias{plotGraph}
\alias{plotColGraph}
\alias{plotRowGraph}
\alias{plotColGraph,ANY,SummarizedExperiment-method}
\alias{plotColGraph,SummarizedExperiment,missing-method}
\alias{plotRowGraph,ANY,SummarizedExperiment-method}
\alias{plotRowGraph,SummarizedExperiment,missing-method}
\title{Plotting igraph objects with information from a \code{SummarizedExperiment}}
\usage{
plotColGraph(x, y, ...)

plotRowGraph(x, y, ...)

\S4method{plotColGraph}{ANY,SummarizedExperiment}(
  x,
  y,
  show_label = FALSE,
  add_legend = TRUE,
  layout = "kk",
  edge_type = c("fan", "link", "arc", "parallel"),
  edge_colour_by = NULL,
  edge_width_by = NULL,
  colour_by = NULL,
  shape_by = NULL,
  size_by = NULL,
  by_exprs_values = "counts",
  other_fields = list(),
  ...
)

\S4method{plotColGraph}{SummarizedExperiment,missing}(x, y, name = "graph", ...)

\S4method{plotRowGraph}{ANY,SummarizedExperiment}(
  x,
  y,
  show_label = FALSE,
  add_legend = TRUE,
  layout = "kk",
  edge_type = c("fan", "link", "arc", "parallel"),
  edge_colour_by = NULL,
  edge_width_by = NULL,
  colour_by = NULL,
  shape_by = NULL,
  size_by = NULL,
  by_exprs_values = "counts",
  other_fields = list(),
  ...
)

\S4method{plotRowGraph}{SummarizedExperiment,missing}(x, y, name = "graph", ...)
}
\arguments{
\item{x, y}{a graph object and a
\code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
object or just a
\code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}.
For the latter object a graph object must be stored in \code{metadata(x)$name}.}

\item{...}{additional arguments for plotting. See
\code{\link{mia-plot-args}} for more details i.e. call \code{help("mia-plot-args")}}

\item{show_label}{\code{logical} (scalar), \code{integer} or \code{character}
vector. If a \code{logical} scalar is given, should tip labels be plotted
or if a logical vector is provided, which labels should be shown? If an
\code{integer} or \code{character} vector is provided, it will be converted
to a logical vector. The \code{integer} values must be in the range of 1
and number of nodes, whereas the values of a \code{character} vector must
match values of a \code{label} or \code{name} column in the node data. In
case of a \code{character} vector only values corresponding to actual
labels will be plotted and if no labels are provided no labels will be
shown. (default: \code{show_label = FALSE})}

\item{add_legend}{logical scalar. Should legends be plotted?
(default: \code{add_legend = TRUE})}

\item{layout}{layout for the plotted graph. See
\code{\link[ggraph:ggraph]{ggraph}} for details. (default:
\code{layout = "kk"})}

\item{edge_type}{type of edge plotted on the graph. See
\code{\link[ggraph:geom_edge_fan]{geom_edge_fan}} for details and other
available geoms. (default:
\code{edge_type = "fan"})}

\item{edge_colour_by}{Specification of a edge metadata field to use for
setting colours of the edges.}

\item{edge_width_by}{Specification of a edge metadata field to use for
setting width of the edges.}

\item{colour_by}{Specification of a column metadata field or a feature to
colour graph nodes by, see the by argument in
\code{\link[scater:retrieveCellInfo]{?retrieveCellInfo}} for possible
values.}

\item{shape_by}{Specification of a column metadata field or a feature to
shape graph nodes by, see the by argument in
\code{\link[scater:retrieveCellInfo]{?retrieveCellInfo}} for possible
values.}

\item{size_by}{Specification of a column metadata field or a feature to
size graph nodes by, see the by argument in
\code{\link[scater:retrieveCellInfo]{?retrieveCellInfo}} for possible
values.}

\item{by_exprs_values}{A string or integer scalar specifying which assay to
obtain expression values from, for use in point aesthetics - see the
\code{exprs_values} argument in
\code{\link[scater:retrieveCellInfo]{?retrieveCellInfo}}.}

\item{other_fields}{Additional fields to include in the node information
without plotting them.}

\item{name}{If \code{x} is a
\code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
the key for subsetting the \code{metadata(x)} to a graph object.}
}
\value{
a \code{\link{ggtree}} plot
}
\description{
\code{plotGraph} plots an \code{igraph} object with additional information
matched from a \code{SummarizedExperiment} object for the nodes only.
Information on the edges have to provided manually.
}
\details{
:
Internally \code{tidygraph} and \code{ggraph} are used. Therefore, all
graph types which can be converted by \code{tidygraph::as_tbl_graph} can
be used.
}
\examples{
\donttest{
# data setup
library(mia)
data(GlobalPatterns)
data(col_graph)
data(row_graph)
data(row_graph_order)
metadata(GlobalPatterns)$col_graph <- col_graph

genus <- mergeFeaturesByRank(GlobalPatterns,"Genus",na.rm=TRUE)
metadata(genus)$row_graph <- row_graph
order <- mergeFeaturesByRank(genus,"Order",na.rm=TRUE)
metadata(order)$row_graph <- row_graph_order

# plot a graph independently
plotColGraph(col_graph,
             genus,
             colour_by = "SampleType",
             edge_colour_by = "weight",
             edge_width_by = "weight",
             show_label = TRUE)

# plot the graph stored in the object
plotColGraph(genus,
             name = "col_graph",
             colour_by = "SampleType",
             edge_colour_by = "weight",
             edge_width_by = "weight")
             

# plot a graph independently
plotRowGraph(row_graph,
             genus,
             colour_by = "Kingdom",
             edge_colour_by = "weight",
             edge_width_by = "weight")

# plot the graph stored in the object
plotRowGraph(genus,
             name = "row_graph",
             colour_by = "Phylum",
             edge_colour_by = "weight",
             edge_width_by = "weight")

                           
# plot a graph independently
plotRowGraph(row_graph_order,
             order,
             colour_by = "Kingdom",
             edge_colour_by = "weight",
             edge_width_by = "weight")

# plot the graph stored in the object and include some labels
plotRowGraph(order,
             name = "row_graph",
             colour_by = "Phylum",
             edge_colour_by = "weight",
             edge_width_by = "weight", 
             show_label = c("Sulfolobales","Spirochaetales",
                            "Verrucomicrobiales"))
                            
# labels can also be included via selecting specific rownames of x/y
plotRowGraph(order,
             name = "row_graph",
             colour_by = "Phylum",
             edge_colour_by = "weight",
             edge_width_by = "weight", 
             show_label = c(1,10,50))
             
# labels can also be included via a logical vector, which has the same length
# as nodes are present
label_select <- rep(FALSE,nrow(order))
label_select[c(1,10,50)] <-  TRUE
plotRowGraph(order,
             name = "row_graph",
             colour_by = "Phylum",
             edge_colour_by = "weight",
             edge_width_by = "weight",
             show_label = label_select)
}
}
