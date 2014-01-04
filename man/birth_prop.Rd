\name{birth_prop}
\alias{birth_prop}
\title{Calculate Bx, the proportion of births for each age x and state}
\usage{
birth_prop(birth_rates, transition_probs, survivor_prop)
}
\arguments{
  \item{birth_rates}{a vector of age-specific birth rates,
  stacked per state}

  \item{transition_probs}{list of transtion probability
  matrices Mx, orderd by age x}

  \item{survivor_prop}{list of survivorship proportions
  matrices Sx, orderd by age x}
}
\value{
list of Bx matrices
}
\description{
\code{birth_prop} returns a list with \code{Bx} matrices
containing the proportion of births, based on the birth
rates \code{Fx}, the transition probabilities \code{Px} and
the survivorship proportion \code{Sx}.
}
\keyword{manip}

