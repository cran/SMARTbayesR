#' @title Simultaneous Upper Credible Intervals
#'
#' @description Compute simultaneous upper credible intervals from draws of embedded dynamic treatment regime probabilities given by the function argument ``thetadraws''.
#'
#' @param thetadraws draws of the embedded dynamic treatment regimes.
#' @param alpha the probability of excluding the true best EDTR from the set of best.
#' @param design specifies to which SMART design to apply function: either design-1, general, or design-3.
#' @param type summary statistic: log-OR, log-RR, or RD
#' @return Upper 1-alpha level simultaneous credible interval limits for the embedded dynamic treatment regimes.
#'
#' @examples
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
#'                design = "design-1",
#'                type = "log-OR")
#' @export

MCBUpperLimits <- function(thetadraws,
                           alpha = 0.05,
                           design = "design-1",
                           type="log-OR") {
  
  if (!(design %in% c("design-1", "general", "design-3"))) {
    stop("'design' must be 'design-1', 'general', or 'design-3'")
  }
  if (!(type %in% c("log-OR","log-RR","RD"))) {
    stop("'type' must be 'log-OR', 'log-RR', or 'RD'")
  }
  if (design == "design-1" & ncol(thetadraws) != 4) {
    stop("'thetadraws' must have 4 EDTRs for design-1")
  }
  if (design == "general" & ncol(thetadraws) != 8) {
    stop("'thetadraws' must have 8 EDTRs for general")
  }
  if (design == "design-3" & ncol(thetadraws) != 6) {
    stop("'thetadraws' must have 6 EDTRs for design-3")
  }
  if (alpha <= 0 | alpha >= 0.5) {
    stop("Probability of exclusion from set of best the optimal EDTR must be between 0 and 0.5")
  }
  
  if (design == "design-1") {
    return(MCBUpperLimitsDesign1(thetadraws,alpha,type))
  }

  if (design == "general") {
    return(MCBUpperLimitsGeneral(thetadraws,alpha,type))
  }
  
  if (design == "design-3") {
    return(MCBUpperLimitsDesign3(thetadraws,alpha,type))
  }
}
