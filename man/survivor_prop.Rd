\name{survivor_prop}
\alias{survivor_prop}
\title{Calculate Sx, the proportion of survivors for each age x and state}
\usage{
survivor_prop(years_lived)
}
\arguments{
  \item{years_lived}{list of years lived matrices
  \code{Lx}, ordered by age x}
}
\value{
list of Sx matrices
}
\description{
\code{survivor_prop} returns a list with \code{Sx}
matrices, containing the proportion of survivors, based on
the years lived in each state \code{Lx}.
}
\keyword{manip}

