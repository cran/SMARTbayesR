#' @title Internal Functions
#'
#' @description This script contains functions which are 
#' used by other functions in SMARTbayesR.
#'
#' @keywords internal
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

###
### 1:
PosteriorEDTRProbsDesign3 <- function(x) {
  cbind(x[,"p_1"]*(x[,"s1"])+x[,"p_2"]*(1-(x[,"s1"])),
        x[,"p_1"]*(x[,"s1"])+x[,"p_3"]*(1-(x[,"s1"])),
        x[,"p_4"]*(x[,"s2"])+x[,"p_5"]*(1-(x[,"s2"])),
        x[,"p_4"]*(x[,"s2"])+x[,"p_6"]*(1-(x[,"s2"])),
        x[,"p_7"]*(x[,"s3"])+x[,"p_8"]*(1-(x[,"s3"])),
        x[,"p_7"]*(x[,"s3"])+x[,"p_9"]*(1-(x[,"s3"])))
}

#####
#### 2:
PosteriorTrtSeqProbDesign1 <- function(niter, 
                                       dat) {

  # Arguments:
  # niter: number of draws from the posteriors of response probabilities
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
PosteriorTrtSeqProbGeneral <- function(niter, 
                                       dat) {


  # Arguments:
  # niter: number of draws from the posteriors of response probabilities
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

  # Posterior probability of response at the end of stage-1 for each of the two stage-1 treatments
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
####
#### 2:
PosteriorTrtSeqProbDesign3 <- function(niter, 
                                       dat) {
  
  #Arguments
  #niter: number of draws from the posteriors of response probabilities
  #dat: dataset
  
  #End of study binary response indicator
  y <- dat$y
  
  #Stage-1 treatment assignment indicator
  a1 <- dat$a1
  
  #Stage-2 treatment assignment indicator
  a2 <- dat$a2
  
  #End of stage-1 binary response indicator
  s <- dat$s
  
  results <- matrix(NA, nrow = niter, ncol = 12)
  
  colnames(results) <- c("p_1","p_2","p_3","p_4","p_5","p_6","p_7","p_8","p_9","s1","s2","s3")
  
  p_1_results <- p_2_results <- p_3_results <- p_4_results <- p_5_results <- p_6_results <- p_7_results <- p_8_results <- p_9_results<- rep(0.5,niter)
  
  
  s1_results <- s2_results <- s3_results <- rep(0.5,niter)
  
  
  #Simulate from each of the six treatment sequences, the probability of response at the end of the trial from the posterior
  p_1_results <- stats::rbeta(niter, shape1 = sum(y[a1==1&s==1])+1, sum(a1==1&s==1)-sum(y[a1==1&s==1])+1)
  p_2_results <- stats::rbeta(niter, shape1 = sum(y[a1==1&s==0&a2==1])+1, sum(a1==1&s==0&a2==1)-sum(y[a1==1&s==0&a2==1])+1)
  p_3_results <- stats::rbeta(niter, shape1 = sum(y[a1==1&s==0&a2==-1])+1, sum(a1==1&s==0&a2==-1)-sum(y[a1==1&s==0&a2==-1])+1)
  
  p_4_results <- stats::rbeta(niter, shape1 = sum(y[a1==2&s==1])+1, sum(a1==2&s==1)-sum(y[a1==2&s==1])+1)
  p_5_results <- stats::rbeta(niter, shape1 = sum(y[a1==2&s==0&a2==1])+1, sum(a1==2&s==0&a2==1)-sum(y[a1==2&s==0&a2==1])+1)
  p_6_results <- stats::rbeta(niter, shape1 = sum(y[a1==2&s==0&a2==-1])+1, sum(a1==2&s==0&a2==-1)-sum(y[a1==2&s==0&a2==-1])+1)
  
  p_7_results <- stats::rbeta(niter, shape1 = sum(y[a1==3&s==1])+1, sum(a1==3&s==1)-sum(y[a1==3&s==1])+1)
  p_8_results <- stats::rbeta(niter, shape1 = sum(y[a1==3&s==0&a2==1])+1, sum(a1==3&s==0&a2==1)-sum(y[a1==3&s==0&a2==1])+1)
  p_9_results <- stats::rbeta(niter, shape1 = sum(y[a1==3&s==0&a2==-1])+1, sum(a1==3&s==0&a2==-1)-sum(y[a1==3&s==0&a2==-1])+1)
  
  #Simulate from each of the two first stage treatments the probability of response at the end of stage-1 from the posterior
  s1_results <- stats::rbeta(niter, sum(s[a1==1])+1, sum(a1==1)-sum(s[a1==1])+1)
  s2_results <- stats::rbeta(niter, sum(s[a1==2])+1, sum(a1==2)-sum(s[a1==2])+1)
  s3_results <- stats::rbeta(niter, sum(s[a1==3])+1, sum(a1==3)-sum(s[a1==3])+1)
  
  
  results[,"p_1"] <- p_1_results
  results[,"p_2"]<-p_2_results
  results[,"p_3"] <- p_3_results
  results[,"p_4"]<-p_4_results
  results[,"p_5"]<-p_5_results
  results[,"p_6"]<-p_6_results
  results[,"p_7"]<-p_7_results
  results[,"p_8"]<-p_8_results
  results[,"p_9"]<-p_9_results
  
  results[,"s1"] <- s1_results
  results[,"s2"] <- s2_results
  results[,"s3"] <- s3_results
  
  return(results)
}

#### 3:
MCBUpperLimitsDesign1 <- function(thetadraws,
                                  alpha = 0.05,
                                  type="log-OR") {

  # Arguments:
  # thetadraws: draws of target parameter
  # alpha: Probability of excluding optimal embedded dynamic treatment regime from the set of best
  # type: log-OR, log-RR, RD
  
  if (type == "log-OR") {
    
    upper_limit <- rep(NA, 4)
    
    #Compute log-odds
    thetadraws_log_odds <- log(thetadraws/(1-thetadraws))
    
    #Compute index of best EDTR
    max_odds_ind <- which.max(colMeans(thetadraws_log_odds))
    
    #Compute log-odds ratios between each EDTR and best
    Log_OR_matrix <- thetadraws_log_odds - matrix(thetadraws_log_odds[,max_odds_ind], nrow = nrow(thetadraws), ncol = 4)
    
    #Rank log-OR
    rank_matrix <- apply(Log_OR_matrix[,-max_odds_ind], 2, rank, ties.method = "random")
    
    #Find max rank across EDTRs
    rank_max <- apply(rank_matrix, 1, max)
    
    #Create sorted log-OR
    new_dat <- apply(Log_OR_matrix[,], 2, sort)
    
    #Compute 100(1-alpha)% upper quantile
    ranks_quantile <- ceiling(stats::quantile(rank_max, 1 - alpha))
    
    #Compute upper limit of credible interval. One for each log-OR which determines the set of best.
    upper_limit <- new_dat[ranks_quantile,]
    
  }
  
  if (type == "log-RR") {
    
    upper_limit <- rep(NA, 4)
    
    #Compute log of draws
    thetadraws_log_prob <- log(thetadraws)
    
    #Compute index of best EDTR
    max_log_prob_ind <- which.max(colMeans(thetadraws_log_prob))
    
    #Compute log-prob ratios between each EDTR and best
    Log_RR_matrix <- thetadraws_log_prob - matrix(thetadraws_log_prob[,max_log_prob_ind], nrow=nrow(thetadraws),ncol = 4)
    
    #Rank log-RR
    rank_matrix <- apply(Log_RR_matrix[,-max_log_prob_ind], 2, rank, ties.method = "random")
    
    #Find max rank across EDTRs
    rank_max <- apply(rank_matrix, 1, max)
    
    #Create sorted log-RR
    new_dat <- apply(Log_RR_matrix[,], 2, sort)
    
    #Compute 100(1-alpha)% upper quantile of max ranks
    ranks_quantile <- ceiling(stats::quantile(rank_max, 1-alpha))
    
    #Compute upper limit of credible interval. One for each log-RR which determines the set of best.
    upper_limit <- new_dat[ranks_quantile, ]
  }
  
  if (type == "RD") {
    
    upper_limit <- rep(NA, 4)
    
    #Compute index of best EDTR
    max_prob_ind <- which.max(colMeans(thetadraws))
    
    #Compute difference between each EDTR and best
    RD_matrix <- thetadraws-matrix(thetadraws[,max_odds_ind], nrow = nrow(thetadraws),ncol = 4)
    
    #Rank RD
    rank_matrix <- apply(RD_matrix[,-max_prob_ind],2,rank,ties.method="random")
    
    #Find max rank across EDTRs
    rank_max <- apply(rank_matrix,1,max)
    
    #Create sorted RD
    new_dat <- apply(RD_matrix[,],2,sort)
    
    #Compute 100(1-alpha)% upper quantile
    ranks_quantile <- ceiling(stats::quantile(rank_max,1-alpha))
    
    #Compute upper limit of credible interval. One for each RD which determines the set of best.
    upper_limit <- new_dat[ranks_quantile,]
  }
  return(upper_limit)
}

#####
#### 3:
MCBUpperLimitsGeneral <- function(thetadraws, 
                                  alpha = 0.05,
                                  type="log-OR") {

  # Arguments:
  # thetadraws: draws of embedded dynamic treatment regime end of study response probabilities
  # alpha: probability of excluding the optimal embedded DTR from the set of best
  # type: log-OR, log-RR, or RD
  if (type == "log-OR") {
    
    upper_limit <- rep(NA,8)
    
    #Compute log-odds
    thetadraws_log_odds <- log(thetadraws/(1-thetadraws))
    
    #Compute index of best EDTR
    max_odds_ind <- which.max(colMeans(thetadraws_log_odds))
    
    #Compute log-odds ratios between each EDTR and best
    Log_OR_matrix <- thetadraws_log_odds - matrix(thetadraws_log_odds[, max_odds_ind], nrow = nrow(thetadraws), ncol = 8)
    
    #Rank log-OR
    rank_matrix <- apply(Log_OR_matrix[,-max_odds_ind], 2, rank, ties.method = 'random')
    
    #Find max rank across EDTRs
    rank_max <- apply(rank_matrix, 1, max)
    
    #Create sorted log-OR
    new_dat <- apply(Log_OR_matrix[,], 2, sort)
    
    #Compute 100(1-alpha)% upper quantile
    ranks_quantile <- ceiling(stats::quantile(rank_max, 1 - alpha))
    
    #Compute upper limit of credible interval. One for each log-OR which determines the set of best.
    upper_limit <- new_dat[ranks_quantile,]
  }
  if (type == "log-RR") {
    
    upper_limit <- rep(NA,8)
    
    #Compute log-probs
    thetadraws_log_prob <- log(thetadraws)
    
    #Compute index of best EDTR
    max_log_prob_ind <- which.max(colMeans(thetadraws_log_prob))
    
    #Compute log-RR between each EDTR and best
    Log_RR_matrix <- thetadraws_log_prob - matrix(thetadraws_log_prob[,max_odds_ind], nrow = nrow(thetadraws), ncol = 8)
    
    #Rank log-RR
    rank_matrix <- apply(Log_RR_matrix[,-max_log_prob_ind], 2, rank, ties.method = 'random')
    
    #Find max rank across EDTRs
    rank_max <- apply(rank_matrix, 1, max)
    
    #Create sorted log-RR
    new_dat <- apply(Log_RR_matrix[,], 2, sort)
    
    #Compute 100(1-alpha)% upper quantile
    ranks_quantile <- ceiling(stats::quantile(rank_max, 1 - alpha))
    
    #Compute upper limit of credible interval. One for each log-RR which determines the set of best.
    upper_limit <- new_dat[ranks_quantile,]
  }
  
  
  if (type == "RD") {
    upper_limit <- rep(NA,8)
    
    
    #Compute index of best EDTR
    max_log_prob_ind <- which.max(colMeans(thetadraws))
    
    #Compute RD between each EDTR and best
    RD_matrix <- thetadraws - matrix(thetadraws[, max_odds_ind], nrow = nrow(thetadraws), ncol = 8)
    
    #Rank RD
    rank_matrix <- apply(RD_matrix[,-max_log_prob_ind], 2, rank, ties.method = 'random')
    
    #Find max rank across EDTRs
    rank_max <- apply(rank_matrix, 1, max)
    
    #Create sorted RD
    new_dat <- apply(RD_matrix[,], 2, sort)
    
    #Compute 100(1-alpha)% upper quantile
    ranks_quantile <- ceiling(stats::quantile(rank_max,1-alpha))
    
    #Compute upper limit of credible interval. One for each RD which determines the set of best.
    upper_limit <- new_dat[ranks_quantile,]
  }


  return(upper_limit)
}
####
#### 3:
MCBUpperLimitsDesign3 <- function(thetadraws,
                                         alpha = 0.05,
                                         type = "log-OR") {
  #Arguments:
  #thetadraws: draws of embedded dynamic treatment regime draws
  #alpha: Probability of excluding the true best EDTR from the set of best
  
  if (type == "log-OR") {
    
    upper_limit <- rep(NA,6)
    
    #Compute log-odds
    thetadraws_log_odds <- log(thetadraws/(1-thetadraws))
    
    #Compute index of best EDTR
    max_odds_ind <- which.max(colMeans(thetadraws_log_odds))
    
    #Compute log-odds ratios between each EDTR and best
    Log_OR_matrix <- thetadraws_log_odds - matrix(thetadraws_log_odds[,max_odds_ind], nrow = nrow(thetadraws), ncol = 6)
    
    #Rank log-OR
    rank_matrix <- apply(Log_OR_matrix[,-max_odds_ind], 2, rank, ties.method="random")
    
    #Find max rank across EDTRs
    rank_max <- apply(rank_matrix, 1, max)
    
    #Create sorted log-OR
    new_dat <- apply(Log_OR_matrix[,], 2, sort)
    
    #Compute 100(1-alpha)% upper quantile
    ranks_quantile <- ceiling(stats::quantile(rank_max, 1 - alpha))
    
    #Compute upper limit of credible interval. One for each log-OR which determines the set of best.
    upper_limit <- new_dat[ranks_quantile,]
  }
  
  if (type == "log-RR") {
    
    upper_limit <- rep(NA,6)
    
    #Compute log-probs
    thetadraws_log_prob<- log(thetadraws)
    
    #Compute index of best EDTR
    max_log_prob_ind <- which.max(colMeans(thetadraws_log_prob))
    
    #Compute log-RR between each EDTR and best
    Log_RR_matrix <- thetadraws_log_prob - matrix(thetadraws_log_prob[, max_log_prob_ind], nrow = nrow(thetadraws), ncol = 6)
    
    #Rank log-RR
    rank_matrix <- apply(Log_RR_matrix[,-max_log_prob_ind], 2, rank, ties.method = "random")
    
    #Find max rank across EDTRs
    rank_max <- apply(rank_matrix, 1, max)
    
    #Create sorted log-RR
    new_dat <- apply(Log_RR_matrix[,], 2, sort)
    
    #Compute 100(1-alpha)% upper quantile
    ranks_quantile <- ceiling(stats::quantile(rank_max, 1 - alpha))
    
    #Compute upper limit of credible interval. One for each log-RR which determines the set of best.
    upper_limit <- new_dat[ranks_quantile, ]
  }
  
  if (type == "RD") {
    
    upper_limit <- rep(NA,6)
    
    #Compute index of best EDTR
    max_prob_ind <- which.max(colMeans(thetadraws))
    
    #Compute RD between each EDTR and best
    RD_matrix <- thetadraws - matrix(thetadraws[, max_odds_ind], nrow = nrow(thetadraws), ncol = 6)
    
    #Rank RD
    rank_matrix <- apply(RD_matrix[,-max_prob_ind], 2, rank, ties.method="random")
    
    #Find max rank across EDTRs
    rank_max <- apply(rank_matrix, 1, max)
    
    #Create sorted RD
    new_dat <- apply(RD_matrix[,], 2, sort)
    
    #Compute 100(1-alpha)% upper quantile
    ranks_quantile <- ceiling(stats::quantile(rank_max, 1 - alpha))
    
    #Compute upper limit of credible interval. One for each RD which determines the set of best.
    upper_limit <- new_dat[ranks_quantile,]
  }
  
  return(upper_limit)
}


#####
#### 4:
PowerBayesianDesign1 <- function(sample_size = 100,
                                 response_prob = c(0.5, 0.9, 0.3, 0.7, 0.5, 0.8),
                                 stage_one_trt_one_response_prob = 0.7,
                                 stage_one_trt_two_response_prob = 0.5,
                                 type = "log-OR",
                                 threshold,
                                 alpha = 0.05) {

  # Arguments:
  # sample_size: sample size
  # response_prob: probability of response at end of trial for each treatment sequence
  # stage_one_trt_one_response_prob: end of stage-1 response probability for stage-1 treatment one
  # stage_one_trt_two_response_prob: end of stage-1 response probability for stage-1 treatment two
  # threshold: threshold for exclusion from set of best.
  # alpha: probability of excluding best EDTR from the set of best

  upper_limit <- array(NA, dim = c(1000, 10, 4))
  
  pb <- utils::txtProgressBar(min=0,
                       max=1000,
                       initial = 0,
                       style=3)
  
  for (i in 1:1000) {

    # Stage-1 treatment indicator
    a1_binom <- stats::rbinom(sample_size, 1, 0.5)

    # Stage-2 randomization indicator given stage-1 treatment was 1
    a21_binom <- stats::rbinom(sample_size, 1, 0.5)
    # Stage-2 randomization indicator given stage-2 treatment was -1
    a22_binom <- stats::rbinom(sample_size, 1, 0.5)

    # Posterior probability of being randomized to 1 vs -1 in stage 1
    a1 <- stats::rbeta(1000, sum(a1_binom) + 1, sample_size - sum(a1_binom) + 1)

    # Posterior probability of being randomized to 1 vs. -1 in stage-2
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


      # Transform draws from treatment sequence response probabilities and stage-1 treatment response probabilities to
      # embedded DTR response probabilities using Robins' G-computation method
      thetadraws <- cbind(
        p_1_results * (s1) + p_2_results * ((1 - s1)),
        p_1_results * (s1) + p_3_results * ((1 - s1)),
        p_4_results * (s2) + p_5_results * ((1 - s2)),
        p_4_results * (s2) + p_6_results * ((1 - s2))
      )
      
      if (type=="log-OR") {
        logORThreshold <- LogOR(response_prob,
                                stage_one_trt_one_response_prob,
                                stage_one_trt_two_response_prob,
                                design = "design-1")
        
        thetadraws_log_odds <- log(thetadraws/(1-thetadraws))
        
        max_odds_ind <- which.max(colMeans(thetadraws_log_odds))
        
        log_OR_matrix <- thetadraws_log_odds-matrix(thetadraws_log_odds[,max_odds_ind],nrow=1000,ncol=4)
        
        
        rank_matrix <- apply(log_OR_matrix[,-max_odds_ind],2,rank,ties.method = 'random')
        
        rank_max <- apply(rank_matrix,1,max)
        
        
        new_dat <- apply(log_OR_matrix[,],2,sort)
        
        ranks_quantile <- ceiling(stats::quantile(rank_max,1-alpha))
        
        upper_limit[i,j,] <- new_dat[ranks_quantile,]
      }
      
      if (type=="log-RR") {
        
        LogRRThreshold <- LogRR(response_prob,
                                stage_one_trt_one_response_prob,
                                stage_one_trt_two_response_prob,
                                design = "design-1")
        
        thetadraws_log_prob <- log(thetadraws)
        
        max_odds_ind <- which.max(colMeans(thetadraws_log_prob))
        
        log_RR_matrix <- thetadraws_log_prob-matrix(thetadraws_log_prob[,max_odds_ind],nrow=1000,ncol=4)
        
        
        rank_matrix <- apply(log_RR_matrix[,-max_odds_ind],2,rank,ties.method = 'random')
        
        rank_max <- apply(rank_matrix,1,max)
        
        
        new_dat <- apply(log_RR_matrix[,],2,sort)
        
        ranks_quantile <- ceiling(stats::quantile(rank_max,1-alpha))
        
        upper_limit[i,j,] <- new_dat[ranks_quantile,]
      }
      
      if (type=="RD") {
        RDThreshold <- RD(response_prob,
                          stage_one_trt_one_response_prob,
                          stage_one_trt_two_response_prob,
                          design = "design-1")
        
        max_odds_ind <- which.max(colMeans(thetadraws))
        
        RD_matrix <- thetadraws-matrix(thetadraws[,max_odds_ind],nrow=1000,ncol=4)
        
        
        rank_matrix <- apply(RD_matrix[,-max_odds_ind],2,rank,ties.method = 'random')
        
        rank_max <- apply(rank_matrix,1,max)
        
        
        new_dat <- apply(RD_matrix[,],2,sort)
        
        ranks_quantile <- ceiling(stats::quantile(rank_max,1-alpha))
        
        upper_limit[i,j,] <- new_dat[ranks_quantile,]
      }
    }
    utils::setTxtProgressBar(pb, i)

  }
  
  close(pb)
  if ( type == "log-OR") {
    rejection_indices <- which(abs(logORThreshold)>threshold)
    
  }
  if ( type == "log-RR") {
    rejection_indices <- which(abs(LogRRThreshold)>threshold)
    
  }
  if ( type == "RD") {
    rejection_indices <- which(abs(RDThreshold)>threshold)
    
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
                                 stage_one_trt_one_response_prob = 0.5,
                                 stage_one_trt_two_response_prob = 0.7,
                                 type = "log-OR",
                                 threshold,
                                 alpha = 0.05) {

  # Arguments:
  # sample_size: total sample size in SMART study
  # response_prob: probability of response for each of eight embedded treatment sequences (not EDTRs)
  # stage_one_trt_one_response_prob: probability of response at the end of stage-1 to treatment a1 = 1
  # stage_one_trt_two_response_prob: probability of response at the end of stage-1 to treatment a1 = 0
  # type: log-OR, log-RR, or RD
  # threshold: threshold for excluding from set of best
  # alpha: probability of excluding best EDTR from the set of best

  upper_limit <- array(NA, dim = c(1000, 10, 8))



  # Stage-1 response indicator (1 if responder, 0 otherwise)
  s <- rep(NA, sample_size)

  pb <- utils::txtProgressBar(min=0,
                       max=1000,
                       initial = 0,
                       style=3)
  for (i in 1:1000) {

    # Stage-1 treatment indicator
    a1 <- stats::rbinom(sample_size, 1, 0.5)
    #s[a1 == 1] <- stats::rbinom(length(which(a1 == 1)), 1, stage_one_trt_two_response_prob)
    #s[a1 == 0] <- stats::rbinom(length(which(a1 == 0)), 1, stage_one_trt_one_response_prob)
    s[a1 == 1] <- stats::rbinom(length(which(a1 == 1)), 1, stage_one_trt_one_response_prob)
    s[a1 == 0] <- stats::rbinom(length(which(a1 == 0)), 1, stage_one_trt_two_response_prob)
    
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


      # Compute embedded DTR end of study response probability draws using Robins' G-computation formula
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


      if (type == "log-OR") {
        logORThreshold <- LogOR(response_prob,
                                stage_one_trt_one_response_prob,
                                stage_one_trt_two_response_prob,
                                design = "general")
        
        thetadraws_log_odds <- log(thetadraws/(1-thetadraws))
        
        max_odds_ind <- which.max(colMeans(thetadraws_log_odds))
        
        log_OR_matrix <- thetadraws_log_odds-matrix(thetadraws_log_odds[,max_odds_ind],nrow=1000,ncol=8)
        
        
        rank_matrix <- apply(log_OR_matrix[,-max_odds_ind],2,rank,ties.method = 'random')
        
        rank_max <- apply(rank_matrix,1,max)
        
        
        new_dat <- apply(log_OR_matrix[,],2,sort)
        
        ranks_quantile <- ceiling(stats::quantile(rank_max,1-alpha))
        
        upper_limit[i,j,] <- new_dat[ranks_quantile,]
      }
      if (type == "log-RR") {
        
        LogRRThreshold <- LogRR(response_prob,
                                stage_one_trt_one_response_prob,
                                stage_one_trt_two_response_prob,
                                design = "general")
        
        thetadraws_log_prob <- log(thetadraws)
        
        max_odds_ind <- which.max(colMeans(thetadraws_log_prob))
        
        log_RR_matrix <- thetadraws_log_prob-matrix(thetadraws_log_prob[,max_odds_ind],nrow=1000,ncol=8)
        
        
        rank_matrix <- apply(log_RR_matrix[,-max_odds_ind],2,rank,ties.method = 'random')
        
        rank_max <- apply(rank_matrix,1,max)
        
        
        new_dat <- apply(log_RR_matrix[,],2,sort)
        
        ranks_quantile <- ceiling(stats::quantile(rank_max,1-alpha))
        
        upper_limit[i,j,] <-new_dat[ranks_quantile,]
      }
      
      if (type == "RD") {
        RDThreshold <- RD(response_prob,
                          stage_one_trt_one_response_prob,
                          stage_one_trt_two_response_prob,
                          design = "general")
        
        max_odds_ind <- which.max(colMeans(thetadraws))
        
        RD_matrix <- thetadraws-matrix(thetadraws[,max_odds_ind],nrow=1000,ncol=8)
        
        
        rank_matrix <- apply(RD_matrix[,-max_odds_ind],2,rank,ties.method = 'random')
        
        rank_max <- apply(rank_matrix,1,max)
        
        
        new_dat <- apply(RD_matrix[,],2,sort)
        
        ranks_quantile <- ceiling(stats::quantile(rank_max,1-alpha))
        
        upper_limit[i,j,] <-new_dat[ranks_quantile,]
      }
    }
    utils::setTxtProgressBar(pb, i)
    
  }
  
  close(pb)
  
  if ( type == "log-OR") {
    rejection_indices <- which(abs(logORThreshold)>threshold)
    
  }
  if ( type == "log-RR") {
    rejection_indices <- which(abs(LogRRThreshold)>threshold)
    
  }
  if ( type == "RD") {
    rejection_indices <- which(abs(RDThreshold)>threshold)
    
  }
  if (length(rejection_indices) == 1) {
    return(mean(apply(upper_limit, 3, function(x) x)[, rejection_indices] < 0))
  } else {
    return(mean(apply(apply(upper_limit, 3, function(x) x)[, rejection_indices] < 0, 1, prod)))
  }
}

PowerBayesianDesign3 <- function(sample_size=100,
                                        response_prob=c(0.5,0.9,0.3,0.7,0.5,0.8,0.5,0.5,0.5),
                                        stage_one_trt_one_response_prob = 0.7,
                                        stage_one_trt_two_response_prob = 0.5,
                                        stage_one_trt_three_response_prob = 0.5,
                                        type = "log-OR",
                                        threshold,
                                        alpha = 0.05)
{
  
  # Arguments
  # sample_size: total SMART study sample size
  # response_prob: probability of response for each of six embedded treatment sequences
  # stage_one_trt_one_response_prob: probability of response to stage-1 treatment given initial treatment one
  # stage_one_trt_two_response_prob: probability of response to stage-1 treatment given initial treatment two
  # stage_one_trt_three_response_prob: probability of response to stage-1 treatment given initial treatment three
  # type: log-OR, log-RR, RD
  # threshold: threshold for exclusion from set of best
  # alpha: probability of excluding best EDTR from set of best
  
  pb <- utils::txtProgressBar(min=0,
                       max=1000,
                       initial = 0,
                       style=3)
  
  upper_limit <- array(NA, dim=c(1000,10,6))
  
  for (i in 1:1000) {
    
    #Stage-1 treatment indicator
    a1_multinom <-sample(1:3,size=sample_size,replace=TRUE)
    
    #Stage-2 randomization indicator given stage-1 treatment was 1
    a21_binom <-stats::rbinom(sample_size,1,0.5)
    #Stage-2 randomization indicator given stage-1 treatment was 2
    a22_binom <-stats::rbinom(sample_size,1,0.5)
    #Stage-2 randomization indicator given stage-1 treatment was 3
    a23_binom <-stats::rbinom(sample_size,1,0.5)
    
    

    #Posterior probability of being randomized to 1 vs 2 or 3 in stage 1
    a1_dirichlet <- LaplacesDemon::rdirichlet(1000,alpha =c(sum(a1_multinom==1)+1,
                                                            sum(a1_multinom==2)+1,
                                                            sum(a1_multinom==3)+1))

    a1_1 <- a1_dirichlet[,1]
    #Posterior probability of being randomized to 2 vs 1 or 3 in stage 1
    a1_2 <- a1_dirichlet[,2]
    #Posterior probability of being randomized to 3 vs 1 or 2 in stage 1
    a1_3 <- a1_dirichlet[,3]

    #Posterior probability of being randomized to 1 vs. 0 in stage-2
    #This is for first stage being 1
    a21 <- stats::rbeta(1000, sum(a21_binom)+1,sample_size-sum(a21_binom)+1)
    #This is for first stage being 2
    a22 <- stats::rbeta(1000, sum(a22_binom)+1,sample_size-sum(a22_binom)+1)
    #This is for first stage being 3
    
    a23 <- stats::rbeta(1000, sum(a23_binom)+1,sample_size-sum(a23_binom)+1)
    
    #Indicator of response to stage 1 treatment 1
    s1_binom <- stats::rbinom(floor(sample_size*mean(a1_multinom==1)),1,stage_one_trt_one_response_prob)
    
    #Posterior probability of response to stage-1 treatment 1
    s1 <- stats::rbeta(1000, sum(s1_binom)+1,sum(a1_multinom==1)-sum(s1_binom)+1)
    
    #Indicator of response to stage-1 treatment 2
    s2_binom <- stats::rbinom(floor(sample_size*mean(a1_multinom==2)),1,stage_one_trt_two_response_prob)
    
    #Posterior probability of response to stage-1 treatment 2
    s2 <- stats::rbeta(1000, sum(s2_binom)+1,sum((a1_multinom==2))-sum(s2_binom)+1)
    
    #Indicator of response to stage-1 treatment 3
    s3_binom <- stats::rbinom(floor(sample_size*mean(a1_multinom==3)),1,stage_one_trt_three_response_prob)
    
    #Posterior probability of response to stage-1 treatment 3
    s3 <- stats::rbeta(1000, sum(s3_binom)+1,sum((a1_multinom==3))-sum(s3_binom)+1)
    
    #Response indicator at end of study for each of the embedded treatment sequences.
    y_1_results <- stats::rbinom(ceiling(mean((sample_size*s1*a1_1))), 1,response_prob[1])
    y_2_results <- stats::rbinom(ceiling(mean((sample_size*a21*(1-s1)*a1_1))),1,response_prob[2])
    y_3_results <- stats::rbinom(ceiling(mean((sample_size*(1-a21)*(1-s1)*a1_1))),1,response_prob[3])
    
    y_4_results <- stats::rbinom(ceiling(mean((sample_size*s2*a1_2))),1,response_prob[4])
    y_5_results <- stats::rbinom(ceiling(mean((sample_size*a22*(1-s2)*a1_2))),1,response_prob[5])
    y_6_results <- stats::rbinom(ceiling(mean((sample_size*(1-a22)*(1-s2)*a1_2))),1,response_prob[6])
    
    y_7_results <- stats::rbinom(ceiling(mean((sample_size*s2*a1_3))),1,response_prob[7])
    y_8_results <- stats::rbinom(ceiling(mean((sample_size*a23*(1-s3)*a1_3))),1,response_prob[8])
    y_9_results <- stats::rbinom(ceiling(mean((sample_size*(1-a23)*(1-s3)*a1_3))),1,response_prob[9])
    
    for (j in 1:10) {
      
      #M=1000, number of MC samples
      #j is number of times drawing the phi's
      #i is the number of datasets
      
      #Posterior probability of response for each of the six embedded treatment sequences
      p_1_results <- stats::rbeta(1000, shape1 = sum(y_1_results)+1,length(y_1_results)-sum(y_1_results)+1)
      p_2_results <- stats::rbeta(1000, shape1 = sum(y_2_results)+1,length(y_2_results)-sum(y_2_results)+1)
      p_3_results <- stats::rbeta(1000, shape1 = sum(y_3_results)+1,length(y_3_results)-sum(y_3_results)+1)
      
      p_4_results <- stats::rbeta(1000, shape1 = sum(y_4_results)+1,length(y_4_results)-sum(y_4_results)+1)
      p_5_results <- stats::rbeta(1000, shape1 = sum(y_5_results)+1,length(y_5_results)-sum(y_5_results)+1)
      p_6_results <- stats::rbeta(1000, shape1 = sum(y_6_results)+1,length(y_6_results)-sum(y_6_results)+1)
      
      
      p_7_results <- stats::rbeta(1000, shape1 = sum(y_7_results)+1,length(y_7_results)-sum(y_7_results)+1)
      p_8_results <- stats::rbeta(1000, shape1 = sum(y_8_results)+1,length(y_8_results)-sum(y_8_results)+1)
      p_9_results <- stats::rbeta(1000, shape1 = sum(y_9_results)+1,length(y_9_results)-sum(y_9_results)+1)
      
      #Transform draws from treatment sequence response probabilities and stage-1 treatment response probabilities to 
      #embedded DTR response probabilities using Robins' G-computation method
      thetadraws <- cbind(p_1_results*(s1)+p_2_results*((1-s1)),
                          p_1_results*(s1)+p_3_results*((1-s1)),
                          p_4_results*(s2)+p_5_results*((1-s2)),
                          p_4_results*(s2)+p_6_results*((1-s2)),
                          p_7_results*(s3)+p_8_results*((1-s3)),
                          p_7_results*(s3)+p_9_results*((1-s3)))
      
      #Compute set of best
      
      if (type=="log-OR") {
        logORThreshold <- LogOR(response_prob,
                                       stage_one_trt_one_response_prob,
                                       stage_one_trt_two_response_prob,
                                       stage_one_trt_three_response_prob,
                                       design = "design-3")
        
        thetadraws_log_odds <- log(thetadraws/(1-thetadraws))
        max_odds_ind <- which.max(colMeans(thetadraws_log_odds))
        
        log_OR_matrix <- thetadraws_log_odds-matrix(thetadraws_log_odds[,max_odds_ind],nrow=1000,ncol=6)
        
        
        rank_matrix <- apply(log_OR_matrix[,-max_odds_ind],2,rank,ties.method = 'random')
        
        rank_max <- apply(rank_matrix,1,max)
        
        
        new_dat <- apply(log_OR_matrix[,],2,sort)
        
        ranks_quantile <- ceiling(stats::quantile(rank_max,1-alpha))
        
        upper_limit[i,j,] <-new_dat[ranks_quantile,]
      }
      
      if (type=="log-RR") {
        LogRRThreshold <- LogRR(response_prob,
                                       stage_one_trt_one_response_prob,
                                       stage_one_trt_two_response_prob,
                                       stage_one_trt_three_response_prob,
                                       design = "design-3")
        
        thetadraws_log_prob <- log(thetadraws)
        max_odds_ind <- which.max(colMeans(thetadraws_log_prob))
        
        log_RR_matrix <- thetadraws_log_prob-matrix(thetadraws_log_prob[,max_odds_ind],nrow=1000,ncol=6)
        
        
        rank_matrix <- apply(log_RR_matrix[,-max_odds_ind],2,rank,ties.method = 'random')
        
        rank_max <- apply(rank_matrix,1,max)
        
        
        new_dat <- apply(log_RR_matrix[,],2,sort)
        
        ranks_quantile <- ceiling(stats::quantile(rank_max,1-alpha))
        
        upper_limit[i,j,] <-new_dat[ranks_quantile,]
      }
      
      if (type=="RD") {
        RDThreshold <- RD(response_prob,
                          stage_one_trt_one_response_prob,
                          stage_one_trt_two_response_prob,
                          stage_one_trt_three_response_prob,
                          design = "design-3")
        
        max_odds_ind <- which.max(colMeans(thetadraws))
        
        RD_matrix <- thetadraws-matrix(thetadraws[,max_odds_ind],nrow=1000,ncol=6)
        
        
        rank_matrix <- apply(RD_matrix[,-max_odds_ind],2,rank,ties.method = 'random')
        
        rank_max <- apply(rank_matrix,1,max)
        
        
        new_dat <- apply(RD_matrix[,],2,sort)
        
        ranks_quantile <- ceiling(stats::quantile(rank_max,1-alpha))
        
        upper_limit[i,j,] <-new_dat[ranks_quantile,]
      }
    }
    

    
    
    utils::setTxtProgressBar(pb, i)
    
  }
  close(pb)
  
  if ( type == "log-OR") {
    rejection_indices <- which(abs(logORThreshold)>threshold)
    
  }
  if ( type == "log-RR") {
    rejection_indices <- which(abs(LogRRThreshold)>threshold)
    
  }
  if ( type == "RD") {
    rejection_indices <- which(abs(RDThreshold)>threshold)
    
  }
  if (length(rejection_indices)==1){
    return(mean(apply(upper_limit, 3, function(x) x)[, rejection_indices] < 0))
    
  } else {
    return(mean(apply(apply(upper_limit, 3, function(x) x)[, rejection_indices] < 0, 1, prod)))
  }
  
}