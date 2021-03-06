% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_helpers.R
\name{nichevol_labels}
\alias{nichevol_labels}
\title{Labels to represent changes of niche characteristics between ancestors and descendants}
\usage{
nichevol_labels(tree, whole_rec_table, ancestor_line = FALSE,
  present = "1", absent = "0", unknown = "?",
  present_col = "#252525", unknown_col = "#d9d9d9",
  no_change_col = "#b2df8a", retraction_col = "#984ea3",
  expansion_col = "#4daf4a", width = 1, height = 1)
}
\arguments{
\item{tree}{an object of class "phylo".}

\item{whole_rec_table}{matrix of reconstructed bins for nodes and species
derived from a process of maximum parsimony or maximum likelihood reconstruction.
See functions \code{\link{bin_par_rec}} or \code{\link{bin_ml_rec}}.}

\item{ancestor_line}{controls whether ancestor line is plotted.
Default = FALSE.}

\item{present}{(character) code indicating environmental bins in which the
species is present. Default = "1".}

\item{absent}{(character) code indicating environmental bins in which the
species is absent. Default = "0".}

\item{unknown}{(character) code indicating environmental bins in which the
species presence is unknown (uncertain). Default = "?".}

\item{present_col}{color for line representing environments where the species
is present. Default = "#252525".}

\item{unknown_col}{color for line representing environments where the species
presence is unknown (uncertain). Default = "#d9d9d9".}

\item{no_change_col}{color for area of the bar representing environments where
no change has been detected. Default = "#b2df8a".}

\item{retraction_col}{color for area of the bar representing environments where
niche retraction has been detected. Default = "#984ea3".}

\item{expansion_col}{color for area of the bar representing environments where
niche expansion has been detected. Default = "#4daf4a".}

\item{width}{value defining the width of bars representing changes in niches;
default = 1.}

\item{height}{value defining the height of bars representing changes in niches;
default = 1.}
}
\description{
nichevol_labels helps in adding bar-type labels that represent how
species' niches changed from ancestors to descendants.
}
\details{
For the moment, only plots of type "phylogram" with "rightwards" or "leftwards"
directions, created with the function \code{\link[ape]{plot.phylo}} from the
package \code{ape} are supported.

Evolution of ecological niches is represented in one environmental dimension,
with vertical bars indicating if the niche of the descendant has expanded,
retracted, or has not changed compared to its ancestor's niche. Lower values of
environmental variables are represented in the lower part of the bar, and the
opposite part of the bar represents higher values.

Changes in niches (evolution) are defined as follows:
\itemize{
\item if (ancestor == present & descendant == absent) {change <- "retraction"}
\item if (ancestor == present & descendant == present) {change <- "no_change"}
\item if (ancestor == present & descendant == unknown) {change <- "no_change"}
\item if (ancestor == absent & descendant == present) {change <- "expansion"}
\item if (ancestor == absent & descendant == absent) {change <- "no_change"}
\item if (ancestor == absent & descendant == unknown) {change <- "no_change"}
\item if (ancestor == unknown & descendant == absent) {change <- "no_change"}
\item if (ancestor == unknown & descendant == present) {change <- "no_change"}
\item if (ancestor == unknown & descendant == unknown) {change <- "no_change"}
}

If \code{ancestor_line} is TRUE, the ancestor line will be plotted on the bar
representing niche evolution. The line will represent where, in the range of
environmental conditions, the ancestor was present, and where its presence is
uncertain (unknown).
}
\examples{
# a simple tree
data("tree5", package = "nichevol")

# a matrix of niche charactes (1 = present, 0 = absent, ? = unknown)
dataTable <- cbind("241" = rep("1", length(tree5$tip.label)),
                   "242" = rep("1", length(tree5$tip.label)),
                   "243" = c("1", "1", "0", "0", "0"),
                   "244" = c("1", "1", "0", "0", "0"),
                   "245" = c("1", "?", "0", "0", "0"))
rownames(dataTable) <- tree5$tip.label

# list with two objects (tree and character table)
treeWdata <- geiger::treedata(tree5, dataTable)

# Maximum parsimony reconstruction
rec_tab <- smooth_rec(bin_par_rec(treeWdata))

# plotting and adding labels
ape::plot.phylo(tree5, label.offset = 0.04)
nichevol_labels(tree5, rec_tab, height = 0.6)
}
