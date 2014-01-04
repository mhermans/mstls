\name{transfer_matrix}
\alias{transfer_matrix}
\title{Construct Mx, the matrix of observed transfer rates}
\usage{
transfer_matrix(observed_rates, multiple = TRUE, absorbing_state = "last")
}
\arguments{
  \item{observed_rates}{dataframe or matrix with the
  observed rates, stacked per state (see Details)}

  \item{multiple}{allow for multiple transitions between
  states during an interval (default TRUE)}

  \item{absorbing_state}{Indicate which state, if any, is
  the absorbing state (currently not used)}
}
\value{
list of Mx matrices
}
\description{
\code{transfer_matrix} returns a list with age-specific
transfer-rate matrices Mx, based on the observed death and
outmigration rates per state.
}
\keyword{manip}

