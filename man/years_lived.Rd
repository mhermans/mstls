\name{years_lived}
\alias{years_lived}
\title{Calculate Lx, number of years lived in each state by unit birth cohort}
\usage{
years_lived(expected_survivors, transfer_rates)
}
\arguments{
  \item{expected_survivors}{list of number of survivor
  matrices \code{lx}, ordered by age x}

  \item{transfer_rates}{list of transfer-rate matrices Mx,
  orderd by age x}
}
\value{
list of Lx matrices
}
\description{
\code{years_lived} returns a list with \code{Lx} matrices,
containing the number of years lived in each state by unit
birth cohort, based on the number of survivors \code{lx}
and the transfer matrices \code{Mx}.
}
\keyword{manip}

