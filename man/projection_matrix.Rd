\name{projection_matrix}
\alias{projection_matrix}
\title{Construct the generalized Leslie matrix G}
\usage{
projection_matrix(birth_prop, survivor_prop)
}
\arguments{
  \item{birth_prop}{list of birth proportions matrices Bx,
  orderd by age x}

  \item{survivor_prop}{list of survivorship proportions
  matrices Sx, orderd by age x}
}
\value{
matrix G
}
\description{
\code{projection_matrix} combines the birth (\code{Bx}) and
survivorship (\code{Sx}) proportions into a generalized
Leslie matrix \code{G}.
}
\keyword{manip}

