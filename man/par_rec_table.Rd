% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_doc.R
\docType{data}
\name{par_rec_table}
\alias{par_rec_table}
\title{Example of table with results from parsimony reconstructions}
\format{
A character matrix with 11 rows and 20 columns.
}
\usage{
par_rec_table
}
\description{
A character table representing species ecological niches derived from
previous preparation processes and reconstructed niches for ancestors.
Each row represents a species or a node and each column a binary
character in which one or more values of the environmental variable
are categorized as used "1", non used "0", or uncertain "?".
}
\examples{
data("par_rec_table", package = "nichevol")

head(par_rec_table)
}
\keyword{datasets}
