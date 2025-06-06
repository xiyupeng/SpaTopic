% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class.R
\docType{class}
\name{SpaTopic-class}
\alias{SpaTopic-class}
\title{A class of the output from 'SpaTopic'}
\description{
Outputs from function \code{\link{SpaTopic_inference}}. 
A \code{\link{list}} contains the following members:
\itemize{
 \item \code{$Perplexity}. The perplexity is for the training data. 
 Let N be the total number of cells across all images.
   \eqn{Perplexity = exp(-loglikelihood/N)}
 \item \code{$Deviance}. \eqn{Deviance = -2loglikelihood}.
 \item \code{$loglikelihood}. The model log-likelihood. 
 \item \code{$loglike.trace}. The log-likelihood for every collected posterior sample. 
 NULL if \code{trace} = FALSE.
 \item \code{$DIC}. Deviance Information Criterion. NULL if \code{trace} = FALSE.
 \item \code{$Beta}. Topic content matrix with rows as celltypes and columns as topics
 \item \code{$Theta}. Topic prevalent matrix with rows as regions and columns as topics
 \item \code{$Ndk}. Number of cells per topic (col) per region (row).
 \item \code{$Nwk}. Number of cells per topic (col) per celltype (row).
 \item \code{$Z.trace}. Number of times cell being assigned to each topic across all posterior samples.
 We can further compute the posterior distributions of \code{Z} (topic assignment) for 
 individual cells.
 \item \code{$doc.trace}. \code{Ndk} for every collected posterior sample.
  NULL if \code{trace} = FALSE.
 \item \code{$word.trace}. \code{Nwk} for every collected posterior sample.
  NULL if \code{trace} = FALSE.
 \item \code{$cell_topics}. Final topic assignments Z for individual cells.
 \item \code{$parameters}. Model parameters used in the analysis.
}
}
\seealso{
\code{\link{SpaTopic_inference}}
}
