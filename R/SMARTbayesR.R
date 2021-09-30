#' SMARTbayesR: A package for Bayesian computation of optimal embedded dynamic treatment regimes and sample size determination with binary outcomes
#'
#' The SMARTbayesR package allows computation of a set of optimal embedded dynamic treatment regimes for a SMART. 
#' Furthermore, it allows power to be calculated for sample size determination.
#'
#' @section SMARTbayesR functions:
#'
#' SimDesign1 simulates a design-1 type SMART
#' 
#' PosteriorTrtSeqProb draws from the posterior of the probabilities of response for each embedded treatment sequence and stage-1 response probabilities.
#'
#' PosteriorEDTRProbs converts treatment sequence end of study response probabilities, stage-1 response probabilities into end of study embedded dynamic treatment regime response probabilities.
#'
#' MCBUpperLimits calculates simultaneous credible intervals which determines the set of optimal dynamic treatment regimes
#'
#' LogOR computes the log-OR between each embedded dynamic treatment regime and the best.
#'
#' LogRR computes the log-RR between each embedded dynamic treatment regime and the best.
#' 
#' RD computes the risk difference between each embedded dynamic treatment regime and the best.
#'
#' PowerBayesian computes the power in a SMART to exclude embedded dynamic treatment regimes inferior to the best by a specified amount.
#'
#' Please see Artman (2020) <arXiv:2008.02341> for details about the methodology.
#'
#' @docType package
#' @name SMARTbayesR
NULL
# > NULL
