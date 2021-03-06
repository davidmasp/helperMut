% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/indels.R
\name{del_microhomology_VR}
\alias{del_microhomology_VR}
\title{Get microhomology of deletions}
\usage{
del_microhomology_VR(vr, genome = genome_selector())
}
\arguments{
\item{vr}{a VRanges object with deletions}

\item{genome}{a BSgenome object with tha assembly}
}
\value{
a VRanges object with added meta columns
}
\description{
Get microhomology of deletions
}
