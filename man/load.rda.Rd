% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TCGA_functions.R
\name{load.rda}
\alias{load.rda}
\title{Load data from an rda file}
\usage{
load.rda(path.to.rda)
}
\arguments{
\item{path.to.rda}{path to rda file}
}
\value{
TCGA.raw.counts
}
\description{
load rda file that is saved from a download.TCGA.RNAseq function call
}
\examples{
load.rda("/home/user/Desktop/TCGA.raw.counts.rda")
}
