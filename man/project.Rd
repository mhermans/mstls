\name{project}
\alias{project}
\title{Project an intial multistate population forward using G}
\usage{
project(init, pmat, nsteps, lbls = NULL)
}
\arguments{
  \item{init}{initial population vector, stacked per state}

  \item{pmat}{projection matrix/generalized Leslie matrix
  G}

  \item{nsteps}{number of steps to project forward}

  \item{lbls}{optional labels for the age groups}
}
\value{
matrix with the population vectors in the columns
}
\description{
\code{project} takes a intitial multistate population
vector and projects it forward using the generalized Leslie
matrix \code{G}.
}
\keyword{manip}

