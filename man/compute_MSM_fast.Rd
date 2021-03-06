% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/muts.R
\name{compute_MSM_fast}
\alias{compute_MSM_fast}
\title{Compute mutational subtype matrix}
\usage{
compute_MSM_fast(
  vr,
  k = 1,
  sep = ">",
  genome = genome_selector(),
  simplify_set = c("C", "A")
)
}
\arguments{
\item{vr}{A VRanges object with SNV calls.}

\item{k}{The extension of the mutation subtype. This is equivalent to the positions upstream and downstream to be included, k = 1 for trinucleotydes, k = 2 for pentanucleotides.}

\item{sep}{A string to separate the mutation context to the alternative base.}

\item{genome}{A BSgenome object installed in the local computer.}

\item{simplify_set}{Set of base pairs used to simplify the mutation calls.}
}
\value{
A count matrix with samples in rows and mutation subtypes in columns
}
\description{
It takes a SNVs VR object and it outputs the matrix with the mutation
counts for each mutation subtype. Samples are in rows, subtypes in
in columns.

see also \link[SomaticSignatures]{motifMatrix}
}
