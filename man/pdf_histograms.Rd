% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{pdf_histograms}
\alias{pdf_histograms}
\title{Helper function to create PDF files with histograms}
\usage{
pdf_histograms(env_data, occ_data, y_values, sp_names, variable_name,
               CL_lines, limits, col, output_directory)
}
\arguments{
\item{env_data}{list of environmental values in M for all species.}

\item{occ_data}{list of environmental values in occurrences for all species.}

\item{y_values}{list of values for the y axis to be used to represent where
occurrences are distributed across the environmental values in M.}

\item{sp_names}{(character) names of the species for which the process will
be performed.}

\item{variable_name}{(character) name of the variable to be plotted.}

\item{CL_lines}{(numeric) confidence limits to be plotted in the histograms.}

\item{limits}{numeric matrix containing the actual values for the confidence
limits of M.}

\item{col}{color for lines representing the confidence limits of M.}

\item{output_directory}{(character) name of the folder in which results will be
written.}
}
\value{
A PDF file written in the output directory containing all resulting figures.
}
\description{
Helper function to create PDF files with histograms
}
\examples{
# example data
e_data <- list(rnorm(1000, 15, 7), rnorm(800, 20, 6), rnorm(1000, 12, 3))
o_data <- list(sample(seq(5, 29, 0.1), 45), sample(seq(10, 33, 0.1), 40),
               sample(seq(1, 16, 0.1), 50))
for (i in 1:3) {
  names(e_data[[i]]) <- e_data[[i]]
  names(o_data[[i]]) <- o_data[[i]]
}
y_val <- list(rep(3, length(o_data)), rep(4, length(o_data)),
              rep(2, length(o_data)))
s_names <- c("sp1", "sp2", "sp3")
lims <- rbind(c(3.5, 26.47), c(10.83, 29.66), c(6.92, 16.91))

tmpd <- file.path(tempdir(), "Hist_to_check") # temporal directory
dir.create(tmpd)

# the running (before running, create output_directory in current directory)
bins <- pdf_histograms(env_data = e_data, occ_data = o_data, y_values = y_val,
                       sp_names = s_names, variable_name = "Temperature",
                       CL_lines = 95, limits = lims, col = "green",
                       output_directory = tmpd)
}
