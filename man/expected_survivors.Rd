\name{expected_survivors}
\alias{expected_survivors}
\title{Calculate lx, the expected number of survivors at exact age x}
\usage{
expected_survivors(transition_probs, radix = 1e+05)
}
\arguments{
  \item{transition_probs}{list of transition probabilities
  Px, ordered by age x}

  \item{radix}{radix (default 100000)}
}
\value{
list of lx matrices
}
\description{
\code{expected_survivors} returns a list with \code{lx}
matrices, containing the expected number of survivors at
exact age x, based on the list of transition probability
matrices \code{Px}.
}
\keyword{manip}

