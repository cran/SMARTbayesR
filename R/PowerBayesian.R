#' @title Power Computation when Planning a SMART
#'
#' @description This function computes the power for a sequential multiple assignment randomized trial (SMART) of one of two designs: "design-1" or "general".
#'
#' @param sample_size the total SMART study sample size.
#' @param response_prob a vector of probabilities of response for each of embedded treatment sequences.
#'  In the case of design 1, there are 6 and for general design there are 8.
#' @param stage_one_trt_one_response_prob the probability of response to stage-1 treatment for first initial treatment.
#' @param stage_one_trt_two_response_prob the probability of response to stage-1 treatment for second initial treatment.
#' @param rejection_indices a vector of indices of which embedded dynamic treatment regimes to exclude in power calculation.
#' @param design specifies for which SMART design to calculate the power: design-1 or general.
#' @param alpha type I error rate (probability of excluding optimal embedded dynamic treatment regime)
#' @return The power to exclude embedded dynamic treatment regimes specified by rejection_indices from the set of best.
#'
#' @examples
#'
#' \donttest{
#' PowerBayesian(
#'   design = "design-1",
#'   sample_size = 100,
#'   response_prob = c(0.5, 0.9, 0.3, 0.7, 0.5, 0.8),
#'   stage_one_trt_one_response_prob = 0.7,
#'   stage_one_trt_two_response_prob = 0.5,
#'   rejection_indices = 1:2
#' )
#'
#' PowerBayesian(
#'   design = "general",
#'   sample_size = 250,
#'   response_prob = c(0.5, 0.9, 0.7, 0.2, 0.3, 0.8, 0.4, 0.7),
#'   stage_one_trt_one_response_prob = 0.7,
#'   stage_one_trt_two_response_prob = 0.5,
#'   rejection_indices = c(1, 2, 4, 5, 6)
#' )
#' }
#' @export


PowerBayesian <- function(design = "design-1",
                          sample_size = 100,
                          response_prob = c(0.5, 0.9, 0.3, 0.7, 0.5, 0.8),
                          stage_one_trt_one_response_prob = 0.7,
                          stage_one_trt_two_response_prob = 0.5,
                          rejection_indices = 2:3,
                          alpha = 0.05) {
  if (design == "design-1") {
    return(PowerBayesianDesign1(
      sample_size,
      response_prob,
      stage_one_trt_one_response_prob,
      stage_one_trt_two_response_prob,
      rejection_indices,
      alpha
    ))
  }

  if (design == "general") {
    return(PowerBayesianGeneral(
      sample_size,
      response_prob,
      stage_one_trt_one_response_prob,
      stage_one_trt_two_response_prob,
      rejection_indices,
      alpha
    ))
  }
}
