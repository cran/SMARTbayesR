#' @title Risk Differences for Embedded Dynamic Treatment Regimes
#'
#' @description Computes the embedded dynamic treatment regime specific risk differences with the best.
#'
#' @param response_prob the probability of response for each of embedded treatment sequences. 
#' In the case of the design 1 SMART, there are 6, for the general design there are 8, and for design-3, there are 9.
#' @param stage_one_trt_one_response_prob the probability of response to stage-1 treatment given stage-1 treatment one.
#' @param stage_one_trt_two_response_prob the probability of response to stage-1 treatment given stage-1 treatment two.
#' @param stage_one_trt_three_response_prob the probability of response to stage-1 treatment given stage-1 treatment three (for design-3 only)
#' @param design which SMART design: design-1 or general.

#' @return The embedded dynamic treatment regime specific risk differences.
#'
#' @export


RD <- function(response_prob = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8),
               stage_one_trt_one_response_prob = 0.6,
               stage_one_trt_two_response_prob = 0.3,
               stage_one_trt_three_response_prob =0.4,
               design = "general") {
  
  if (any(response_prob < 0 | response_prob > 1)) {
    stop("Treatment sequence response probabilities must be between 0 and 1")
  }
  
  if (stage_one_trt_one_response_prob < 0 | stage_one_trt_one_response_prob > 1 | stage_one_trt_two_response_prob < 0 | stage_one_trt_two_response_prob > 1 | stage_one_trt_three_response_prob < 0 | stage_one_trt_three_response_prob > 1) {
    stop("Stage-1 treatment response probabilities must be between 0 and 1.")
  }
  
  if (!(design %in% c("design-1", "general", "design-3"))) {
    stop("'design' must be 'design-1', 'general', or 'design-3'")
  }
  
  if (design == "design-1" & length(response_prob) != 6) {
    stop("Design-1 must have 6 treatment sequences.")
  }
  
  if (design == "general" & length(response_prob) != 8) {
    stop("The general SMART must have 8 treatment sequences.")
  }
  
  if (design == "design-3" & length(response_prob) != 9) {
    stop("The design-3 SMART must have 9 treatment sequences.")
  }
  
  if (design == "design-1") {
    EDTRs <- c(
      response_prob[1] * stage_one_trt_one_response_prob + response_prob[2] * (1 - stage_one_trt_one_response_prob),
      response_prob[1] * stage_one_trt_one_response_prob + response_prob[3] * (1 - stage_one_trt_one_response_prob),
      response_prob[4] * stage_one_trt_two_response_prob + response_prob[5] * (1 - stage_one_trt_two_response_prob),
      response_prob[4] * stage_one_trt_two_response_prob + response_prob[6] * (1 - stage_one_trt_two_response_prob)
    )
  }
  
  if (design == "general") {
    EDTRs <- c(
      response_prob[1] * stage_one_trt_one_response_prob + response_prob[3] * (1 - stage_one_trt_one_response_prob),
      response_prob[1] * stage_one_trt_one_response_prob + response_prob[4] * (1 - stage_one_trt_one_response_prob),
      response_prob[2] * stage_one_trt_one_response_prob + response_prob[3] * (1 - stage_one_trt_one_response_prob),
      response_prob[2] * stage_one_trt_one_response_prob + response_prob[4] * (1 - stage_one_trt_one_response_prob),
      response_prob[5] * stage_one_trt_two_response_prob + response_prob[7] * (1 - stage_one_trt_two_response_prob),
      response_prob[5] * stage_one_trt_two_response_prob + response_prob[8] * (1 - stage_one_trt_two_response_prob),
      response_prob[6] * stage_one_trt_two_response_prob + response_prob[7] * (1 - stage_one_trt_two_response_prob),
      response_prob[6] * stage_one_trt_two_response_prob + response_prob[8] * (1 - stage_one_trt_two_response_prob)
    )
  }
  
  if (design == "design-3") {
    EDTRs <- c(
      response_prob[1] * stage_one_trt_one_response_prob + response_prob[2] * (1 - stage_one_trt_one_response_prob),
      response_prob[1] * stage_one_trt_one_response_prob + response_prob[3] * (1 - stage_one_trt_one_response_prob),
      response_prob[4] * stage_one_trt_two_response_prob + response_prob[5] * (1 - stage_one_trt_two_response_prob),
      response_prob[4] * stage_one_trt_two_response_prob + response_prob[6] * (1 - stage_one_trt_two_response_prob),
      response_prob[7] * stage_one_trt_three_response_prob + response_prob[8] * (1 - stage_one_trt_three_response_prob),
      response_prob[7] * stage_one_trt_three_response_prob + response_prob[9] * (1 - stage_one_trt_three_response_prob)
    )
  }
  

  ## Compute Risk difference
  thetadraws_RD <- (EDTRs)
  # Compute index of best EDTR
  max_RD_ind <- (which.max((thetadraws_RD)))
  
  # Compute risk differences between each EDTR and best
  RD_output1 <- matrix((thetadraws_RD - thetadraws_RD[max_RD_ind]), nrow = 1, ncol = length(thetadraws_RD))
  
  
  
  
  if (design == "general") {
    colnames(RD_output1) <- c(
      "EDTR 1",
      "EDTR 2",
      "EDTR 3",
      "EDTR 4",
      "EDTR 5",
      "EDTR 6",
      "EDTR 7",
      "EDTR 8"
    )
  }
  
  if (design == "design-1") {
    colnames(RD_output1) <- c(
      "EDTR 1",
      "EDTR 2",
      "EDTR 3",
      "EDTR 4"
    )
  }
  
  if (design == "design-3") {
    colnames(RD_output1) <- c(
      "EDTR 1",
      "EDTR 2",
      "EDTR 3",
      "EDTR 4",
      "EDTR 5",
      "EDTR 6"
    )
  }
  
  return(RD_output1)
}
