#' @title Log-Odds Ratios for Embedded Dynamic Treatment Regimes
#'
#' @description Computes the embedded dynamic treatment regime specific log odds ratios.
#'
#' @param response_prob the probability of response for each of embedded treatment sequences. In the case of the design 1 SMART, there are 6 and for the general design there are 8.
#' @param stage_one_trt_one_response_prob the probability of response to stage-1 treatment given initial treatment one.
#' @param stage_one_trt_two_response_prob the probability of response to stage-1 treatment given initial treatment two.
#' @param design which SMART design: design-1 or general.

#' @return The embedded dynamic treatment regime specific log-OR.
#'
#' @export


LogOR <- function(response_prob = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8),
                  stage_one_trt_one_response_prob = 0.6,
                  stage_one_trt_two_response_prob = 0.3,
                  design = "general") {
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


  # Compute log-OR
  thetadraws_log_odds <- log(EDTRs / (1 - EDTRs))

  # Compute index of best EDTR
  max_odds_ind <- (which.max((thetadraws_log_odds)))

  Log_OR_output1 <- matrix((thetadraws_log_odds - thetadraws_log_odds[max_odds_ind]), nrow = 1, ncol = length(thetadraws_log_odds))

  if (design == "general") {
    colnames(Log_OR_output1) <- c(
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
    colnames(Log_OR_output1) <- c(
      "EDTR 1",
      "EDTR 2",
      "EDTR 3",
      "EDTR 4"
    )
  }
  return(data.table::data.table(Log_OR_output1))
}
