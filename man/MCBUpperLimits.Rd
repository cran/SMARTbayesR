% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MCBUpperLimits.R
\name{MCBUpperLimits}
\alias{MCBUpperLimits}
\title{Simultaneous Upper Credible Intervals}
\usage{
MCBUpperLimits(thetadraws, alpha = 0.05, design = "design-1", type = "log-OR")
}
\arguments{
\item{thetadraws}{draws of the embedded dynamic treatment regimes.}

\item{alpha}{the probability of excluding the true best EDTR from the set of best.}

\item{design}{specifies to which SMART design to apply function: either design-1, general, or design-3.}

\item{type}{summary statistic: log-OR, log-RR, or RD}
}
\value{
Upper 1-alpha level simultaneous credible interval limits for the embedded dynamic treatment regimes.
}
\description{
Compute simultaneous upper credible intervals from draws of embedded dynamic treatment regime probabilities given by the function argument ``thetadraws''.
}
\examples{
dat <- SimDesign1(sample_size=250,
                  response_prob = c(0.5,0.9,0.3,0.7,0.5,0.8),
                  stage_one_trt_one_response_prob = 0.7,
                  stage_one_trt_two_response_prob = 0.4)

x <- PosteriorTrtSeqProb(niter = 1000, dat, design = "design-1")

thetadraws <- PosteriorEDTRProbs(x, design = "design-1")

MCBUpperLimits(thetadraws,
               alpha = 0.05,
               design = "design-1",
               type = "log-OR")
}
