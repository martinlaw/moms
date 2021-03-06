% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findDesNewCP.R
\name{interimDecision}
\alias{interimDecision}
\title{Obtain interim decision for multi-outcome, two-stage drop-the-loser trial}
\usage{
interimDecision(findDTL.output, test.statistics, return.lookup = FALSE)
}
\arguments{
\item{findDTL.output}{The output from a call to findDTL.}

\item{test.statistics}{A vector of observed interim test statistics}

\item{return.lookup}{Logical. Will return a lookup table if TRUE. Default is FALSE.}
}
\value{
Returns a list containing at least two elements. The first element, decision, is a statement regarding whether
the trial should proceed, and if so, what outcomes should be retained for the second stage. The second
element, cp, is a vector of the conditional power values for each outcome.
If return.lookup==TRUE, the function will return a third list element, lookup, which is a lookup table containing the test statistics for each outcome that correspond to a range of CP values.
}
\description{
This function gives the interim decision for multi-outcome, two-stage drop-the-loser designs that declare trial success
when a specified number of outcomes show promise.
}
\examples{
dtl.out <- findDTL(K=4, Kmax=3, m=2, vars=c(1, 1.01, 2, 1.5), delta0=0.1, delta1=0.4, alpha.k=0.05, cp.l=0.3, cp.u=0.95, n.min=10, n.max=40, power=0.8, corr.scalar=0.4, working.outs=c(1,2))
interimDecision(dtl.out, c(0.32, -1.20, 2.01, 1.45))
}
