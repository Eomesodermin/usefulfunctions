% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TCGA_functions.R
\name{preprocess.TCGA}
\alias{preprocess.TCGA}
\title{Preprocess TCGA RNAseq data}
\usage{
preprocess.TCGA(
  input.data,
  dataset.name = NULL,
  quantile.filt = TRUE,
  save.normed.data = TRUE,
  output.dir = "data/"
)
}
\arguments{
\item{input.data}{output from download.TCGA.RNAseq function}

\item{dataset.name}{character string containing the name of the dataset}

\item{quantile.filt}{logical indicating whether to filter the data by quantile}

\item{save.normed.data}{logical indicating whether to save the normalised data to file}

\item{output.dir}{character string containing the directory to save the normalised data to}
}
\value{
data.frame containing the normalised data
}
\description{
Process the output of download.TCGA.RNAseq function by filtering for protein coding genes only, normalising by GC content, and quantile filtering for low expressed genes
}
\examples{
raw.data <- download.TCGA.RNAseq("TCGA-BRCA")
or
raw.data <- load.rda("path.to.rda.file")

normed.data <- preprocess.TCGA(input.data = raw.data, dataset.name = "TCGA_BRCA")

}
