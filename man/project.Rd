\name{project}
\alias{project}
\title{Project an intial multistate population forward using G}
\usage{
project(init_pop, n_states, pmat, n_steps, lbls = NULL)
}
\arguments{
  \item{init_pop}{initial population vector, stacked per
  state}

  \item{n_states}{number of states}

  \item{pmat}{projection matrix/generalized Leslie matrix
  G}

  \item{n_steps}{number of steps to project forward}

  \item{lbls}{optional labels for the age groups}
}
\value{
matrix with the stacked population vectors in the columns
}
\description{
\code{project} takes a intitial multistate population
vector and projects it forward using the generalized Leslie
matrix \code{G}.
}
\keyword{manip}

