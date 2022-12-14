% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/me2_BS_read_correlations.R
\name{me2_BS_read_correlations}
\alias{me2_BS_read_correlations}
\title{Get the factor "correlations" from the BS results file}
\usage{
me2_BS_read_correlations(
  me2_bs_txt_file,
  tidy_output = FALSE,
  block_boundaries = list(start = "factor \\"correlations\\" with Best-fit factors")
)
}
\arguments{
\item{me2_bs_txt_file}{file name of the BS results file.}

\item{tidy_output}{Should the output be reshaped into tidy data? Default:
FALSE}

\item{block_boundaries}{A list containing the \dQuote{start} string used to
identify the boundaries of the block containing the values for the
correlations. Typically, this block is at the end of the file and therefor
no \dQuote{end} string is required.}
}
\value{
tibble containing the factor "correlations" with Best-fit factors
}
\description{
Typical use of Bootstrap using Multilinear Engine version 2 (ME-2) will
provide three different files after a run: \dQuote{.dat}: machine readable
format, \dQuote{.rsd}: results for the residuals and \dQuote{.txt}: text file
with auxiliary information (i.e. headers). The headers in the text (.txt)
file are used in me2tools to split the data into several blocks. In this
function the factor "correlations" with the Best-fit factors are read using
the provided header information, denoting the \dQuote{start} line of the
block containing the required data.
}
