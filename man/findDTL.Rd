% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findDesNewCP.R
\name{findDTL}
\alias{findDTL}
\title{Find Multi-Outcome Two-Stage Drop-the-Loser designs}
\usage{
findDTL(
  K,
  Kmax,
  m,
  alpha.k,
  power,
  corr.mat = NULL,
  vars = NULL,
  corr.scalar = NULL,
  delta0,
  delta1,
  delta.true = NULL,
  cp.l,
  cp.u,
  n.min,
  n.max,
  working.outs = NULL,
  nsims = 1000
)
}
\arguments{
\item{K}{Number of outcomes}

\item{Kmax}{Maximum number of outcomes permitted in stage 2}

\item{m}{Number of outcomes required to show promise for trial success}

\item{alpha.k}{The desired type-I error-rate.}

\item{power}{The desired power.}

\item{corr.mat}{A square matrix of the correlations between outcomes. Must be K-dimensional and have 1's on the diagonal.}

\item{vars}{A vector of outcome variances. If single value is entered, it is used for all outcomes with a warning.}

\item{corr.scalar}{A scalar of the correlation between outcomes. If entered, it is used for all correlations with a warning.}

\item{delta0}{A vector of anticipated lower effect sizes. If a single value is entered, it is used for all outcomes with a warning.}

\item{delta1}{A vector of anticipated upper effect sizes. If a single value is entered, it is used for all outcomes with a warning.}

\item{delta.true}{Optional. A matrix of true effect sizes (with number of columns==K). If only 2 columns are supplied, will take delta.true\link{1} as true delta for all working outcomes and delta.true\link{2} as true delta for all non-working outcomes.}

\item{cp.l}{The lower bound for conditional power.}

\item{cp.u}{The upper bound for conditional power.}

\item{n.min}{The minimum sample size to search over.}

\item{n.max}{The maximum sample size to search over.}

\item{working.outs}{A vector of the indices of outcomes that are taken to be the "working" or "best-performing" outcomes for the purposes of calculating the sample size. If not given, the first m outcomes will be used, with a warning.}

\item{nsims}{The number of trials simulated. Default is 1000.}
}
\value{
The function returns a list of length two The first element, input, contains the values inputted into the call.
The second element, results, gives the final and interim stopping boundaries and the operating characteristics.
}
\description{
This function finds multi-outcome, two-stage drop-the-loser designs that declare trial success
when a specified number of outcomes show promise. This function uses simulation.
}
\details{
if delta.true is used, an additional list element is returned, true.results, containing the operating characteristics of the obtained design taking into account the true effect sizes supplied.
}
\examples{
findDTL(K=4, Kmax=3, m=2, vars=c(1, 1.01, 2, 1.5), delta0=0.1, delta1=0.4, alpha.k=0.05, cp.l=0.3, cp.u=0.95, n.min=10, n.max=40, power=0.8, corr.scalar=0.4, working.outs=c(1,2))

m1 <- matrix(NA, 4, 4)
m1[lower.tri(m1, diag=F)] <- vec
m1 <- t(m1)
m1[lower.tri(m1, diag=F)] <- vec
diag(m1) <- 1
findDTL(nsims = 1e3, K=4, Kmax=3, m=2, vars = c(1, 1.01, 2, 1.5), delta0 = 0.1, delta1 = 0.4, alpha.k = 0.05, cp.l = 0.3, cp.u = 0.95, n.min = 10, n.max = 40, power = 0.8, corr.mat = m1, working.outs=c(1,2))
}
