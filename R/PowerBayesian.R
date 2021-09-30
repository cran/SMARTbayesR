#' @title Power Calculation for a SMART with a Binary Outcome
#'
#' @description This function computes the power for a sequential multiple assignment randomized trial (SMART) of one of three designs: "design-1" or "general" or "design-3".
#'
#' @param sample_size the total SMART study sample size.
#' @param response_prob a vector of probabilities of response for each of embedded treatment sequences.
#'  In the case of design 1, there are 6, for general design there are 8, and for design-3 there are 9
#' @param stage_one_trt_one_response_prob the probability of response to stage-1 treatment for first stage-1 treatment.
#' @param stage_one_trt_two_response_prob the probability of response to stage-1 treatment for second stage-1 treatment.
#' @param stage_one_trt_three_response_prob the probability of response to stage-1 treatment for third stage-1 treatment (for design-3 only).
#' @param design specifies for which SMART design to calculate the power: design-1, general, or design-3.
#' @param type specifies log-OR, RD or log-RR.
#' @param threshold minimum detectable difference between each EDTR and the best
#' @param alpha probability of excluding optimal embedded dynamic treatment regime
#' @return The power to exclude embedded dynamic treatment regimes bigger than threshold from the set of best.
#'
#' @examples
#' \donttest{
#' PowerBayesian(
#'   design = "design-1",
#'   sample_size = 100,
#'   response_prob = c(0.5, 0.9, 0.3, 0.7, 0.5, 0.8),
#'   stage_one_trt_one_response_prob = 0.7,
#'   stage_one_trt_two_response_prob = 0.5,
#'   type="log-OR",
#'   threshold=0.2
#' )
#'
#' PowerBayesian(
#'   design = "general",
#'   sample_size = 250,
#'   response_prob = c(0.5, 0.9, 0.7, 0.2, 0.3, 0.8, 0.4, 0.7),
#'   stage_one_trt_one_response_prob = 0.7,
#'   stage_one_trt_two_response_prob = 0.5,
#'   type="log-OR",
#'   threshold=0.2
#' )
#' }
#' @export


PowerBayesian <- function(design = "design-1",
                          sample_size = 100,
                          response_prob = c(0.5, 0.9, 0.3, 0.7, 0.5, 0.8),
                          stage_one_trt_one_response_prob = 0.7,
                          stage_one_trt_two_response_prob = 0.5,
                          stage_one_trt_three_response_prob = 0.4,
                          type = "log-OR",
                          threshold, 
                          alpha = 0.05) {
  if (threshold <= 0) {
    stop("Threshold must be positive")
  }
  if (sample_size <= 0 | sample_size != round(sample_size)) {
    stop("Sample size n must be a positive integer")
  }
  if (alpha <=0 | alpha >= 0.5) {
    stop("Probability of exclusion from set of best the optimal EDTR must be between 0 and 0.5")
  }
  if (!(type %in% c("log-OR","log-RR","RD"))) {
    stop("'type' must be 'log-OR', 'log-RR', or 'RD'")
  }
  if (any(response_prob < 0|response_prob > 1)) {
    stop("Treatment sequence response probabilities must be between 0 and 1")
  }
  
  if (stage_one_trt_one_response_prob < 0 | stage_one_trt_one_response_prob > 1 | stage_one_trt_two_response_prob < 0 | stage_one_trt_two_response_prob > 1 | stage_one_trt_three_response_prob < 0 | stage_one_trt_three_response_prob > 1) {
    stop("Stage-1 treatment response probabilities must be between 0 and 1")
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
    return(PowerBayesianDesign1(
      sample_size,
      response_prob,
      stage_one_trt_one_response_prob,
      stage_one_trt_two_response_prob,
      type,
      threshold, 
      alpha
    ))
  }

  if (design == "general") {
    return(PowerBayesianGeneral(
      sample_size,
      response_prob,
      stage_one_trt_one_response_prob,
      stage_one_trt_two_response_prob,
      type,
      threshold, 
      alpha
    ))
  }
  
  if (design == "design-3") {
    return(PowerBayesianDesign3(
      sample_size,
      response_prob,
      stage_one_trt_one_response_prob,
      stage_one_trt_two_response_prob,
      stage_one_trt_three_response_prob,
      type,
      threshold, 
      alpha
    ))
  }
}
