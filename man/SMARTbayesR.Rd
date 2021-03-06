% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SMARTbayesR.R
\docType{package}
\name{SMARTbayesR}
\alias{SMARTbayesR}
\title{SMARTbayesR: A package for Bayesian computation of optimal dynamic treatment regimes and sample size determination}
\description{
The SMARTbayesR package allows computation of a set of optimal embedded dynamic treatment regimes for a SMART. Furthermore, it allows power to be calculated for sample size determination.
}
\details{
Five of the important functions are as follows:
PosteriorTrtSeqProb, PosteriorEDTRProbs, MCBUpperLimits, LogOR, and PowerBayesian.
}
\section{SMARTbayesR functions}{


PosteriorTrtSeqProb draws from the posterior of the probabilities of response for each treatment sequence and stage-1 response probabilities.

PosteriorEDTRProbs converts treatment sequence end of study response probabilities, stage-1 response probabilities into end of study DTR response probabilities.

MCBUpperLimits calculates simultaneous credible intervals which determines the set of optimal DTRs.

LogOR computes the log-OR between each embedded dynamic treatment regime and the best.

PowerBayesian computes the power in a SMART to exlude DTRs inferior to the best by a specified amount.

Please see Artman (2020) <arXiv:2008.02341> for details about the methodology.
}

