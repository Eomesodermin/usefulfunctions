% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Gene_Ontology.R
\name{GO.function}
\alias{GO.function}
\title{Gene ontology}
\usage{
GO.function(
  markers,
  topn = 1000,
  org = "human",
  title.var = "",
  plot = TRUE,
  ontology.type = "BP",
  ...
)
}
\arguments{
\item{markers}{a dataframe with a column "gene" and a column "FDR"}

\item{topn}{number of top genes to be used for GO enrichment}

\item{org}{character for which organism the genes come from either "human" or "mouse"}

\item{title.var}{character to use as title of the plot}

\item{plot}{logical, if TRUE, plot the result, if FALSE, return the result}

\item{ontology.type}{string for which ontology type to assess one of either "BP", "MF", or "CC"}

\item{...}{other parameters to be passed to enrichGO function}
}
\value{
if plot is TRUE, return a plot, if plot is FALSE, return a dataframe
}
\description{
Function is designed to take only up or downregulated genes
Function will find GO terms enriched in the list of up or downregulated genes
input requires a column "gene" and a column "FDR"
gene should be as SYMBOL
}
\author{
Dillon Corvino wrote the wrapper for clusterProfiler function
}
