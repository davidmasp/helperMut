% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/profiles.R
\name{cos_sim_vect}
\alias{cos_sim_vect}
\title{Cosine similarity between 2 named vectors}
\usage{
cos_sim_vect(x, y)
}
\arguments{
\item{x}{A vector (optionally named)}

\item{y}{A vector (optionally named)}
}
\value{
A single double indicating the cosine distant
}
\description{
It reorder the given vectors according to the names in x and performs a
cosine similarity computation.
}
\examples{

cos_sim_vect(x = c(1,2,3),c(2,3,4))

}
\references{
See here \url{https://en.wikipedia.org/wiki/Cosine_similarity}
}
\seealso{
Similar function based on matrices available at \link[lsa]{cosine}
}
