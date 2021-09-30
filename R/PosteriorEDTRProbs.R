#' @title Convert Treatment Sequence Draws into Embedded Dynamic Treatment Regime Draws
#'
#' @description Apply Robins' G-computation formula to compute the embedded dynamic treatment regime draws from as a weighted average of treatment sequence and stage-1 response probability draws.
#'
#' If design is "design-1", then compute for Design 1 SMART with 6 embedded treatment sequences and 4 embedded dynamic treatment regimes.
#'
#' If design is "general", then compute for General SMART with 8 embedded treatment sequences and 8 embedded dynamic treatment regimes.
#'
#' If design is "design-3", then compute for Design 3 SMART with 9 embedded treatment sequences and 6 embedded dynamic treatment regimes.
#' 
#' @param x A data frame consisting of draws from the posterior of the end of study response probabilities of each treatment sequence 
#' and of stage-1 response probabilities for each stage-1 treatment
#' @param design Which SMART design to compute the posterior draws for: "design-1" or "general" or "design-3".
#' @return Matrix of EDTR specific posterior response probability draws at the end of the study
#' There will be 4 columns for design-1, 8 columns for design general, and 6 columns for design-3 each corresponding to an EDTR. The number of rows will be the same as that of x.
#'
#' @examples
#' dat <- SimDesign1(sample_size=250,
#'                   response_prob = c(0.5,0.9,0.3,0.7,0.5,0.8),
#'                   stage_one_trt_one_response_prob = 0.7,
#'                   stage_one_trt_two_response_prob = 0.4)
#'
#' x <- PosteriorTrtSeqProb(niter = 1000, dat, design = "design-1")
#'
#' PosteriorEDTRProbs(x, design = "design-1")
#' @export
#' @details
#' For the General SMART design, x should have columns
#' p_1,
#' p_2,
#' p_3,
#' p_4,
#' p_5,
#' p_6,
#' p_7,
#' p_8,
#' s1, and s2.
#'
#'
#'
#' For the Design-1 SMART, x should have columns
#' p_1,
#' p_2,
#' p_3,
#' p_4,
#' p_5,
#' p_6,
#' s1, and s2.
#' 
#' For the design-3 SMART, x should have columns
#' 
#' p_1,
#' p_2,
#' p_3,
#' p_4,
#' p_5,
#' p_6,
#' p_7,
#' p_8,
#' p_9
#' s1, s2, and s3
#' 
#' These are the posterior draws of the response probabilities for each treatment sequence and stage-1 response probability draws.
#'
#' s1 contains the draws of the stage-1 response probability for the first treatment, s2 is analogous for the second treatment and s3 for the third treatment.

PosteriorEDTRProbs <- function(x, design = "design-1") {
  
  
  if (!(design %in% c("design-1", "general", "design-3"))) {
    stop("'design' must be 'design-1', 'general', or 'design-3'")
  }
  
  if (design == "design-1") {
    return(PosteriorEDTRProbsDesign1(x))
  }

  if (design == "general") {
    return(PosteriorEDTRProbsGeneral(x))
  }
  
  if (design == "design-3") {
    return(PosteriorEDTRProbsDesign3(x))
  }
}
