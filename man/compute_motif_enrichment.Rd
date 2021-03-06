% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{compute_motif_enrichment}
\alias{compute_motif_enrichment}
\title{Compute mutation motif enrichemnt}
\usage{
compute_motif_enrichment(
  vr,
  genome,
  k_offset = 20,
  set = "TCW>K",
  control_set = "NCN>D"
)
}
\arguments{
\item{vr}{a VRanges object with mutation data}

\item{genome}{a BSgenome object}

\item{k_offset}{k for window offset, used to control for ntp at risk}

\item{set}{a mutation set to test}

\item{control_set}{a mutation set to control for}
}
\value{
A dataframe with tidied results of the fisher test. The dataframe contains
the enrichment value as a column.
}
\description{
Given a VRanges object which can contain either a sample, a set of samples
or a fraction of a sample, calculate the enrichment of a motif compared
to a control motif.
}
\details{
A given window size (k_offset*2+1) is used to compute ntps at risk.

Implemets the the enrichment measure used in
Roberts Nature Genetics (2013) for apobec enrichment in mutation
clusters.

In brief, the enrichment is calculated as:

\deqn{E = ( muts_set * ctx_control ) / (muts_control * ctx_set)}

Note that this value is similar to the odds ratio, however, in the
fisher test and its odds ratio, the control set excludes the test set
while it's included in the enrichment measure.
}
