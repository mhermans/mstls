\name{transfer_prob}
\alias{transfer_prob}
\title{Calculate age-specific state transfer probabilities Px.}
\usage{
transfer_prob(transfer_rates, multiple = TRUE)
}
\arguments{
  \item{transfer_rates}{list of transfer matrices Mx,
  orderd by age x}

  \item{multiple}{allow for multiple transitions between
  states during an interval (default TRUE)}
}
\value{
list of Px matrices
}
\description{
\code{transfer_prob} returns a list of age-specific
transfer probabilities matrices Px, based on observed
transfer rates (death and outmigration) structured in the
transition matrices Mx. It allows for single or multiple
transitions per interval.
}
\keyword{manip}

