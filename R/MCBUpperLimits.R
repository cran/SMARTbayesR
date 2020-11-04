#' @title Simultaneous Upper Credible Intervals
#'
#' @description Compute simultaneous upper credible intervals from draws of embedded dynamic treatment regime probabilities given by the function argument ``thetadraws''.
#'
#' @param thetadraws draws of the embedded dynamic treatment regimes.
#' @param alpha the type I error rate (probability of excluding the true best EDTR).
#' @param design specifies to which SMART design to apply function: either design-1 or general.
#'
#' @return Upper 1-alpha level simultaneous credible interval limits for each of the embedded dynamic treatment regimes.
#'
#' @examples
#'
#'
#'dat <- SimDesign1(sample_size=250,
#'                   response_prob = c(0.5,0.9,0.3,0.7,0.5,0.8),
#'                   stage_one_trt_one_response_prob = 0.7,
#'                   stage_one_trt_two_response_prob = 0.4)
#'
#' x <- PosteriorTrtSeqProb(niter = 1000, dat, design = "design-1")
#'
#' thetadraws <- PosteriorEDTRProbs(x, design = "design-1")
#'
#' MCBUpperLimits(thetadraws,
#'                alpha = 0.05,
#'                design = "design-1")
#'
#'
#'
#' @export

MCBUpperLimits <- function(thetadraws, alpha = 0.05, design = "design-1") {
  if (design == "design-1") {
    return(MCBUpperLimitsDesign1(thetadraws,alpha))
  }

  if (design == "general") {
    return(MCBUpperLimitsGeneral(thetadraws,alpha))
  }
}
