#' @title Internal Functions
#'
#' @description This script contains functions which are implicitly
#' used by other functions in SMARTBayesR.
#'
#'
#' @noRd
#'


#####
#### 1:
PosteriorEDTRProbsDesign1 <- function(x) {
  return(cbind(
    x[, "p_1"] * (x[, "s1"]) + x[, "p_2"] * (1 - (x[, "s1"])),
    x[, "p_1"] * (x[, "s1"]) + x[, "p_3"] * (1 - (x[, "s1"])),
    x[, "p_4"] * (x[, "s2"]) + x[, "p_5"] * (1 - (x[, "s2"])),
    x[, "p_4"] * (x[, "s2"]) + x[, "p_6"] * (1 - (x[, "s2"]))
  ))
}

#####
#### 1:
PosteriorEDTRProbsGeneral <- function(x) {
  return(cbind(
    x[, "p_1"] * x[, "s1"] + x[, "p_3"] * (1 - x[, "s1"]),
    x[, "p_1"] * x[, "s1"] + x[, "p_4"] * (1 - x[, "s1"]),
    x[, "p_2"] * x[, "s1"] + x[, "p_3"] * (1 - x[, "s1"]),
    x[, "p_2"] * x[, "s1"] + x[, "p_4"] * (1 - x[, "s1"]),
    x[, "p_5"] * x[, "s2"] + x[, "p_7"] * (1 - x[, "s2"]),
    x[, "p_5"] * x[, "s2"] + x[, "p_8"] * (1 - x[, "s2"]),
    x[, "p_6"] * x[, "s2"] + x[, "p_7"] * (1 - x[, "s2"]),
    x[, "p_6"] * x[, "s2"] + x[, "p_8"] * (1 - x[, "s2"])
  ))
}

#####
#### 2:
PosteriorTrtSeqProbDesign1 <- function(niter, dat) {

  # Arguments:
  # niter: number of iterations
  # dat: dataset

  # End of study binary response indicator
  y <- dat$y

  # Stage-1 treatment assignment indicator
  a1 <- dat$a1

  # Stage-2 treatment assignment indicator
  a2 <- dat$a2

  # End of stage-1 binary response indicator
  s <- dat$s

  results <- matrix(NA, nrow = niter, ncol = 8)

  colnames(results) <- c("p_1", "p_2", "p_3", "p_4", "p_5", "p_6", "s1", "s2")

  p_1_results <- p_2_results <- p_3_results <- p_4_results <- p_5_results <- p_6_results <- rep(0.5, niter)


  s1_results <- s2_results <- rep(0.5, niter)


  # Simulate from each of the six treatment sequences, the probability of response at the end of the trial from the posterior
  p_1_results <- stats::rbeta(niter, shape1 = sum(y[a1 == 1 & s == 1]) + 1, sum(a1 == 1 & s == 1) - sum(y[a1 == 1 & s == 1]) + 1)
  p_2_results <- stats::rbeta(niter, shape1 = sum(y[a1 == 1 & s == 0 & a2 == 1]) + 1, sum(a1 == 1 & s == 0 & a2 == 1) - sum(y[a1 == 1 & s == 0 & a2 == 1]) + 1)
  p_3_results <- stats::rbeta(niter, shape1 = sum(y[a1 == 1 & s == 0 & a2 == -1]) + 1, sum(a1 == 1 & s == 0 & a2 == -1) - sum(y[a1 == 1 & s == 0 & a2 == -1]) + 1)
  p_4_results <- stats::rbeta(niter, shape1 = sum(y[a1 == -1 & s == 1]) + 1, sum(a1 == -1 & s == 1) - sum(y[a1 == -1 & s == 1]) + 1)
  p_5_results <- stats::rbeta(niter, shape1 = sum(y[a1 == -1 & s == 0 & a2 == 1]) + 1, sum(a1 == -1 & s == 0 & a2 == 1) - sum(y[a1 == -1 & s == 0 & a2 == 1]) + 1)
  p_6_results <- stats::rbeta(niter, shape1 = sum(y[a1 == -1 & s == 0 & a2 == -1]) + 1, sum(a1 == -1 & s == 0 & a2 == -1) - sum(y[a1 == -1 & s == 0 & a2 == -1]) + 1)

  # Simulate from each of the two first stage treatments the probability of response at the end of stage-1 from the posterior
  s1_results <- stats::rbeta(niter, sum(s[a1 == 1]) + 1, sum(a1 == 1) - sum(s[a1 == 1]) + 1)
  s2_results <- stats::rbeta(niter, sum(s[a1 == -1]) + 1, sum(a1 == -1) - sum(s[a1 == -1]) + 1)


  results[, "p_1"] <- p_1_results
  results[, "p_2"] <- p_2_results
  results[, "p_3"] <- p_3_results
  results[, "p_4"] <- p_4_results
  results[, "p_5"] <- p_5_results
  results[, "p_6"] <- p_6_results
  results[, "s1"] <- s1_results
  results[, "s2"] <- s2_results

  return(results)
}

#####
#### 2:
PosteriorTrtSeqProbGeneral <- function(niter, dat) {


  # Arguments:
  # niter: number of iterations
  # dat: dataset

  results <- matrix(NA, nrow = niter, ncol = 10)

  colnames(results) <- c("p_1", "p_2", "p_3", "p_4", "p_5", "p_6", "p_7", "p_8", "s1", "s2")

  # End of study binary outcome
  y <- dat$y

  # Stage-1 treatment assignment
  a1 <- dat$a1

  # Stage-2 treatment assignment for responders to stage-1 treatment
  a2r <- dat$a2r

  # Stage-2 treatment assignment for non-responders to stage-1 treatment
  a2nr <- dat$a2nr

  # End of stage-1 response indicator (1 if responder at end of stage-1, 0 otherwise)
  s <- dat$s


  # Posterior probabilities of response at the end of the study for each of the eight embedded treatment sequences
  p_1_results <- stats::rbeta(niter, shape1 = sum(y[a1 == 1 & s == 1 & a2r == 1]) + 1, shape2 = sum(a1 == 1 & s == 1 & a2r == 1) - sum(y[a1 == 1 & s == 1 & a2r == 1]) + 1)
  p_2_results <- stats::rbeta(niter, shape1 = sum(y[a1 == 1 & s == 1 & a2r == -1]) + 1, shape2 = sum(a1 == 1 & s == 1 & a2r == -1) - sum(y[a1 == 1 & s == 1 & a2r == -1]) + 1)
  p_3_results <- stats::rbeta(niter, shape1 = sum(y[a1 == 1 & s == 0 & a2nr == 1]) + 1, shape2 = sum(a1 == 1 & s == 0 & a2nr == 1) - sum(y[a1 == 1 & s == 0 & a2nr == 1]) + 1)
  p_4_results <- stats::rbeta(niter, shape1 = sum(y[a1 == 1 & s == 0 & a2nr == -1]) + 1, shape2 = sum(a1 == 1 & s == 0 & a2nr == -1) - sum(y[a1 == 1 & s == 0 & a2nr == -1]) + 1)
  p_5_results <- stats::rbeta(niter, shape1 = sum(y[a1 == -1 & s == 1 & a2r == 1]) + 1, shape2 = sum(a1 == -1 & s == 1 & a2r == 1) - sum(y[a1 == -1 & s == 1 & a2r == 1]) + 1)
  p_6_results <- stats::rbeta(niter, shape1 = sum(y[a1 == -1 & s == 1 & a2r == -1]) + 1, shape2 = sum(a1 == -1 & s == 1 & a2r == -1) - sum(y[a1 == -1 & s == 1 & a2r == -1]) + 1)
  p_7_results <- stats::rbeta(niter, shape1 = sum(y[a1 == -1 & s == 0 & a2nr == 1]) + 1, shape2 = sum(a1 == -1 & s == 0 & a2nr == 1) - sum(y[a1 == -1 & s == 0 & a2nr == 1]) + 1)
  p_8_results <- stats::rbeta(niter, shape1 = sum(y[a1 == -1 & s == 0 & a2nr == -1]) + 1, shape2 = sum(a1 == -1 & s == 0 & a2nr == -1) - sum(y[a1 == -1 & s == 0 & a2nr == -1]) + 1)

  #  Posterior probability of response at the end of stage-1 for each of the two stage-1 treatments
  s1_results <- stats::rbeta(niter, sum(s[a1 == 1] == 1) + 1, length(which(a1 == 1)) - sum(s[a1 == 1]) + 1)
  s2_results <- stats::rbeta(niter, sum(s[a1 == -1] == 1) + 1, length(which(a1 == -1)) - sum(s[a1 == -1]) + 1)


  results[, "p_1"] <- p_1_results
  results[, "p_2"] <- p_2_results
  results[, "p_3"] <- p_3_results
  results[, "p_4"] <- p_4_results
  results[, "p_5"] <- p_5_results
  results[, "p_6"] <- p_6_results
  results[, "p_7"] <- p_7_results
  results[, "p_8"] <- p_8_results

  results[, "s1"] <- s1_results
  results[, "s2"] <- s2_results
  return(results)
}

#### 3:
MCBUpperLimitsDesign1 <- function(thetadraws,alpha = 0.05) {

  # Arguments:
  # thetadraws: draws of target parameter
  # alpha: type I error rate (excluding optimal embedded dynamic treatment regime)

  upper_limit <- rep(NA, 4)

  # Compute log-OR
  thetadraws_log_odds <- log(thetadraws / (1 - thetadraws))

  # Compute index of best EDTR
  max_odds_ind <- which.max(colMeans(thetadraws_log_odds))

  # Compute log-odds ratios between each EDTR and best
  Log_OR_matrix <- thetadraws_log_odds - matrix(thetadraws_log_odds[, max_odds_ind], nrow = nrow(thetadraws), ncol = 4)

  # Rank log-OR
  rank_matrix <- apply(Log_OR_matrix, 2, rank, ties.method = "min")

  # Find max rank
  rank_max <- apply(rank_matrix, 1, max)

  # Create sorted log-OR
  new_dat <- apply(Log_OR_matrix[rank_max, ], 2, sort)

  # Compute 100(1-alpha)% upper quantile
  ranks_quantile <- ceiling(stats::quantile(rank_max, 1 - alpha))

  # Compute upper limit of credible interval. One for each log-OR which determines the set of best.
  upper_limit <- new_dat[ranks_quantile, ]


  return(upper_limit)
}

#####
#### 3:
MCBUpperLimitsGeneral <- function(thetadraws, alpha = 0.05) {

  # Arguments:
  # thetadraws: draws of embedded dynamic treatment regime end of study response probabilities
  # alpha: type I error rate (excluding the optimal embedded DTR)

  upper_limit <- rep(NA, 8)

  # Compute log-OR
  thetadraws_log_odds <- log(thetadraws / (1 - thetadraws))

  # Compute index of best EDTR
  max_odds_ind <- which.max(colMeans(thetadraws_log_odds))

  # Compute log-odds ratios between each EDTR and best
  Log_OR_matrix <- thetadraws_log_odds - matrix(thetadraws_log_odds[, max_odds_ind], nrow = nrow(thetadraws), ncol = 8)

  # Rank log-OR
  rank_matrix <- apply(Log_OR_matrix, 2, rank, ties.method = "min")

  # Find max rank
  rank_max <- apply(rank_matrix, 1, max)

  # Create sorted log-OR
  new_dat <- apply(Log_OR_matrix[rank_max, ], 2, sort)

  # Compute 100(1-alpha)% upper quantile
  ranks_quantile <- ceiling(stats::quantile(rank_max, 1 - alpha))

  # Compute upper limit of credible interval. One for each log-OR which determines the set of best.
  upper_limit <- new_dat[ranks_quantile, ]


  return(upper_limit)
}


#####
#### 4:
PowerBayesianDesign1 <- function(sample_size = 100,
                                 response_prob = c(0.5, 0.9, 0.3, 0.7, 0.5, 0.8),
                                 stage_one_trt_one_response_prob = 0.7,
                                 stage_one_trt_two_response_prob = 0.5,
                                 rejection_indices = 2:3,
                                 alpha = 0.05) {

  # Arguments:
  # sample_size: sample size
  # response_prob: probability of response at end of trial for each treatment sequence
  # stage_one_trt_one_response_prob: end of stage-1 response probability for stage-1 treatment one
  # stage_one_trt_two_response_prob: end of stage-1 response probability for stage-1 treatment two
  # rejection_indices: which EDTRs to exclude
  # alpha: type I error rate

  upper_limit <- array(NA, dim = c(500, 10, 4))

  for (i in 1:500) {

    # Stage-1 treatment indicator
    a1_binom <- stats::rbinom(sample_size, 1, 0.5)

    # Stage-2 randomization indicator given stage-1 treatment was 1
    a21_binom <- stats::rbinom(sample_size, 1, 0.5)
    # Stage-2 randomization indicator given stage-2 treatment was -1
    a22_binom <- stats::rbinom(sample_size, 1, 0.5)

    # Posterior probability of being randomized to 1 vs 0 in stage 1
    a1 <- stats::rbeta(1000, sum(a1_binom) + 1, sample_size - sum(a1_binom) + 1)

    # Posterior probability of being randomized to 1 vs. 0 in stage-2
    # This is for first stage being 1
    a21 <- stats::rbeta(1000, sum(a21_binom) + 1, sample_size - sum(a21_binom) + 1)
    # This is for first stage being 0
    a22 <- stats::rbeta(1000, sum(a22_binom) + 1, sample_size - sum(a22_binom) + 1)

    # Indicator of response to stage 1 treatment 1
    s1_binom <- stats::rbinom(floor(sample_size * mean(a1_binom)), 1, stage_one_trt_one_response_prob)
    # Posterior probability of response to stage-1 treatment 1
    s1 <- stats::rbeta(1000, sum(s1_binom) + 1, sum(a1_binom) - sum(s1_binom) + 1)

    # Indicator of response to stage-1 treatment 0
    s2_binom <- stats::rbinom(floor(sample_size * mean(1 - a1_binom)), 1, stage_one_trt_two_response_prob)

    # Posterior probability of response to stage-1 treatment 0
    s2 <- stats::rbeta(1000, sum(s2_binom) + 1, sum((1 - a1_binom)) - sum(s2_binom) + 1)

    # Response indicator at end of study for each of the embedded treatment sequences.
    y_1_results <- stats::rbinom(ceiling(mean((sample_size * s1 * a1))), 1, response_prob[1])
    y_2_results <- stats::rbinom(ceiling(mean((sample_size * a21 * (1 - s1) * a1))), 1, response_prob[2])
    y_3_results <- stats::rbinom(ceiling(mean((sample_size * (1 - a21) * (1 - s1) * a1))), 1, response_prob[3])

    y_4_results <- stats::rbinom(ceiling(mean((sample_size * s2 * (1 - a1)))), 1, response_prob[4])
    y_5_results <- stats::rbinom(ceiling(mean((sample_size * a22 * (1 - s2) * (1 - a1)))), 1, response_prob[5])
    y_6_results <- stats::rbinom(ceiling(mean((sample_size * (1 - a22) * (1 - s2) * (1 - a1)))), 1, response_prob[6])

    for (j in 1:10) {
      # M=1000, number of MC samples
      # j=number of times drawing the phi's
      # i is the number of datasets
      # Posterior probability of response for each of the six embedded treatment sequences
      p_1_results <- stats::rbeta(1000, shape1 = sum(y_1_results) + 1, length(y_1_results) - sum(y_1_results) + 1)
      p_2_results <- stats::rbeta(1000, shape1 = sum(y_2_results) + 1, length(y_2_results) - sum(y_2_results) + 1)
      p_3_results <- stats::rbeta(1000, shape1 = sum(y_3_results) + 1, length(y_3_results) - sum(y_3_results) + 1)

      p_4_results <- stats::rbeta(1000, shape1 = sum(y_4_results) + 1, length(y_4_results) - sum(y_4_results) + 1)
      p_5_results <- stats::rbeta(1000, shape1 = sum(y_5_results) + 1, length(y_5_results) - sum(y_5_results) + 1)
      p_6_results <- stats::rbeta(1000, shape1 = sum(y_6_results) + 1, length(y_6_results) - sum(y_6_results) + 1)


      # Transform draws from treatment sequence response probabilites and stage-1 treatment response probabilites to
      # embedded DTR response probabilites using Robin's G-computation method
      thetadraws <- cbind(
        p_1_results * (s1) + p_2_results * ((1 - s1)),
        p_1_results * (s1) + p_3_results * ((1 - s1)),
        p_4_results * (s2) + p_5_results * ((1 - s2)),
        p_4_results * (s2) + p_6_results * ((1 - s2))
      )

      # Apply Betensky's method
      thetadraws_log_odds <- log(thetadraws / (1 - thetadraws))
      max_odds_ind <- which.max(colMeans(thetadraws_log_odds))

      log_OR_matrix <- thetadraws_log_odds - matrix(thetadraws_log_odds[, max_odds_ind], nrow = 1000, ncol = 4)


      rank_matrix <- apply(log_OR_matrix, 2, rank, ties.method = "min")

      rank_max <- apply(rank_matrix, 1, max)


      new_dat <- apply(log_OR_matrix[rank_max, ], 2, sort)

      ranks_quantile <- ceiling(stats::quantile(rank_max, 1 - alpha))

      upper_limit[i, j, ] <- new_dat[ranks_quantile, ]
      # print(i)
    }
  }
  if (length(rejection_indices) == 1) {
    return(mean(apply(upper_limit, 3, function(x) x)[, rejection_indices] < 0))
  } else {
    return(mean(apply(apply(upper_limit, 3, function(x) x)[, rejection_indices] < 0, 1, prod)))
  }
}

#####
#### 4:
# Compute power from relevant inputs
PowerBayesianGeneral <- function(sample_size = 500,
                                 response_prob = c(0.5, 0.9, 0.7, 0.2, 0.3, 0.8, 0.4, 0.7),
                                 stage_one_trt_one_response_prob = 0.7,
                                 stage_one_trt_two_response_prob = 0.5,
                                 rejection_indices = c(1, 2, 4, 5, 6),
                                 alpha = 0.05) {

  # Arguments:
  # sample_size: total sample size in SMART study
  # response_prob: probability of response for each of eight embedded treatment sequences (not EDTRs)
  # stage_one_trt_one_response_prob: probability of response at the end of stage-1 to treatment a1 = 1
  # stage_one_trt_two_response_prob: probability of response at the end of stage-1 to treatment a1 = 0
  # rejection_indices: indices of embedded DTRs to exclude from the set of best in power calculation
  # alpha: Type I error rate

  upper_limit <- array(NA, dim = c(500, 10, 8))



  # Stage-1 response indicator (1 if responder, 0 otherwise)
  s <- rep(NA, sample_size)


  for (i in 1:500) {

    # Stage-1 treatment indicator
    a1 <- stats::rbinom(sample_size, 1, 0.5)
    s[a1 == 1] <- stats::rbinom(length(which(a1 == 1)), 1, stage_one_trt_two_response_prob)
    s[a1 == 0] <- stats::rbinom(length(which(a1 == 0)), 1, stage_one_trt_one_response_prob)

    # Stage-2 randomization indicator for responders to stage-1 treatment
    a2r <- 2 * stats::rbinom(sample_size, 1, 0.5) - 1

    # Stage-2 randomization indicator for non-responders to stage-1 treatment
    a2nr <- 2 * stats::rbinom(sample_size, 1, 0.5) - 1


    # End of study response indicators for each of the eight embedded treatment sequences (not DTRs) (1 is response, 0 is non-response)
    y1 <- stats::rbinom(length(which(a1 == 1 & s == 1 & a2r == 1)), 1, (response_prob[1]))
    y2 <- stats::rbinom(length(which(a1 == 1 & s == 1 & a2r == -1)), 1, (response_prob[2]))
    y3 <- stats::rbinom(length(which(a1 == 1 & s == 0 & a2nr == 1)), 1, (response_prob[3]))
    y4 <- stats::rbinom(length(which(a1 == 1 & s == 0 & a2nr == -1)), 1, (response_prob[4]))
    y5 <- stats::rbinom(length(which(a1 == 0 & s == 1 & a2r == 1)), 1, (response_prob[5]))
    y6 <- stats::rbinom(length(which(a1 == 0 & s == 1 & a2r == -1)), 1, (response_prob[6]))
    y7 <- stats::rbinom(length(which(a1 == 0 & s == 0 & a2nr == 1)), 1, (response_prob[7]))
    y8 <- stats::rbinom(length(which(a1 == 0 & s == 0 & a2nr == -1)), 1, (response_prob[8]))


    for (j in 1:10) {
      # Draw 1000 draws from the posterior of the probability of response at the end of the study for each of the eight embedded treatment sequences.
      p_1_results <- stats::rbeta(1000, shape1 = sum(y1) + 1, shape2 = sum(a1 == 1 & s == 1 & a2r == 1) - sum(y1) + 1)
      p_2_results <- stats::rbeta(1000, shape1 = sum(y2) + 1, shape2 = sum(a1 == 1 & s == 1 & a2r == -1) - sum(y2) + 1)
      p_3_results <- stats::rbeta(1000, shape1 = sum(y3) + 1, shape2 = sum(a1 == 1 & s == 0 & a2nr == 1) - sum(y3) + 1)
      p_4_results <- stats::rbeta(1000, shape1 = sum(y4) + 1, shape2 = sum(a1 == 1 & s == 0 & a2nr == -1) - sum(y4) + 1)
      p_5_results <- stats::rbeta(1000, shape1 = sum(y5) + 1, shape2 = sum(a1 == 0 & s == 1 & a2r == 1) - sum(y5) + 1)
      p_6_results <- stats::rbeta(1000, shape1 = sum(y6) + 1, shape2 = sum(a1 == 0 & s == 1 & a2r == -1) - sum(y6) + 1)
      p_7_results <- stats::rbeta(1000, shape1 = sum(y7) + 1, shape2 = sum(a1 == 0 & s == 0 & a2nr == 1) - sum(y7) + 1)
      p_8_results <- stats::rbeta(1000, shape1 = sum(y8) + 1, shape2 = sum(a1 == 0 & s == 0 & a2nr == -1) - sum(y8) + 1)



      # Draw 1000 draws from the posterior of the probability of response at the end of stage 1 for stage-1 treatment 1 and 0, respectively.
      s1 <- stats::rbeta(1000, sum(s[a1 == 1]) + 1, sum(a1 == 1) - sum(s[a1 == 1]) + 1)
      s2 <- stats::rbeta(1000, sum(s[a1 == 0]) + 1, sum(a1 == 0) - sum(s[a1 == 0]) + 1)


      # Compute embedded DTR end of study response probability draws using Robin's G-computation formula
      thetadraws <- cbind(
        p_1_results * (s1) + p_3_results * (1 - (s1)),
        p_1_results * (s1) + p_4_results * (1 - (s1)),
        p_2_results * (s1) + p_3_results * (1 - (s1)),
        p_2_results * (s1) + p_4_results * (1 - (s1)),
        p_5_results * (s2) + p_7_results * (1 - (s2)),
        p_5_results * (s2) + p_8_results * (1 - (s2)),
        p_6_results * (s2) + p_7_results * (1 - (s2)),
        p_6_results * (s2) + p_8_results * (1 - (s2))
      )


      # Perform Bayesian MCB
      thetadraws_log_odds <- log(thetadraws / (1 - thetadraws))
      max_odds_ind <- which.max(colMeans(thetadraws_log_odds))

      log_OR_matrix <- thetadraws_log_odds - matrix(thetadraws_log_odds[, max_odds_ind], nrow = 1000, ncol = 8)


      rank_matrix <- apply(log_OR_matrix, 2, rank, ties.method = "min")

      rank_max <- apply(rank_matrix, 1, max)


      new_dat <- apply(log_OR_matrix[rank_max, ], 2, sort)

      ranks_quantile <- ceiling(stats::quantile(rank_max, 1 - alpha))

      upper_limit[i, j, ] <- new_dat[ranks_quantile, ]
    }
  }
  if (length(rejection_indices) == 1) {
    return(mean(apply(upper_limit, 3, function(x) x)[, rejection_indices] < 0))
  } else {
    return(mean(apply(apply(upper_limit, 3, function(x) x)[, rejection_indices] < 0, 1, prod)))
  }
}

