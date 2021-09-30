#' @title Treatment Sequence Response Probabilities from Dataset
#'
#' @description Draws from the posterior of the treatment sequence response probabilities.
#'
#' @param niter the number of posterior draws.
#' @param dat a data frame (see Details).
#' @param design an indicator of which SMART design, design-1, general, or design-3.
#' @return Posterior draws of the probability of response at the end of the study for each embedded treatment sequence and
#' the posterior draws of the probability of response at the end of stage-1 for each stage-1 treatment.
#'
#' @examples
#' dat <- SimDesign1(sample_size=250,
#'                               response_prob = c(0.5,0.9,0.3,0.7,0.5,0.8),
#'                               stage_one_trt_one_response_prob = 0.7,
#'                               stage_one_trt_two_response_prob = 0.4)
#'
#' PosteriorTrtSeqProb(niter = 1000, dat, design = "design-1")
#'
#' @export
#' @details
#' dat should contain the following columns:
#'
#' y, the end of study binary response indicator
#'
#' a1, the stage-1 treatment assignment indicator
#'
#' s, the end of stage-1 binary response indicator
#'
#' Additionally, for design-1 and design-3 it should contain
#' a2, the stage-2 treatment assignment indicator
#'
#' For the general design, it should contain
#'
#' a2r, stage-2 treatment assignment for responders to stage-1 treatment.
#'
#' a2nr, stage-2 treatment assignment for non-responders to stage-1 treatment.

PosteriorTrtSeqProb <- function(niter, dat, design = "design-1") {
  
  if (niter <= 0 | niter != round(niter)) {
    stop("Sample size n must be a positive integer")
  }
  
  
  if (!(design %in% c("design-1", "general", "design-3"))) {
    stop("'design' must be 'design-1', 'general', or 'design-3'")
  }
  
  
  if (design == "design-1") {
    return(PosteriorTrtSeqProbDesign1(niter, dat))
  }

  if (design == "general") {
    return(PosteriorTrtSeqProbGeneral(niter, dat))
  }
  
  if (design == "design-3") {
    return(PosteriorTrtSeqProbDesign3(niter, dat))
  }
}
