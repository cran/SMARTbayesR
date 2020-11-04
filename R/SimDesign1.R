#' @title Simulate a SMART with Design 1
#'
#' @description This function simulates a SMART with Design 1.
#'
#' @param sample_size the total sample size
#' @param response_prob a vector of probabilities of response for each of the 6 embedded treatment sequences.
#' @param stage_one_trt_one_response_prob the probability of response to first stage-1 treatment.
#' @param stage_one_trt_two_response_prob the probability of response to second stage-1 treatment.
#'
#' @return a data frame with treatment response indicators for each stage of treatment, a1 and a2,
#' end of stage-1 response indicator s, and final outcome y.
#'
#'
#' @examples
#'
#' dat <- SimDesign1(sample_size=250,
#'                              response_prob = c(0.5,0.9,0.3,0.7,0.5,0.8),
#'                              stage_one_trt_one_response_prob = 0.7,
#'                              stage_one_trt_two_response_prob = 0.4)
#'
#' @export

SimDesign1 <-
  function(sample_size=250,
           response_prob=c(0.5,0.5,0.5,0.8,0.7,0.5),
           stage_one_trt_one_response_prob = 0.7,
           stage_one_trt_two_response_prob = 0.4){





    a1 <- 2*stats::rbinom(sample_size,1,.5)-1 #First stage treatment indicator: coded as -1 and +1

    a2 <- 2*stats::rbinom(sample_size,1,0.5)-1 #Second stage treatment indicator: coded as -1 and +1

    #Stage-1 response probabilities
    s<-rep(NA,sample_size)
    s[a1==1] <- stats::rbinom(length(which(a1==1)),size=1,stage_one_trt_one_response_prob)
    s[a1==-1] <- stats::rbinom(length(which(a1==-1)),size=1,stage_one_trt_two_response_prob)

    #End-of-study outcomes
    y <- rep(NA,sample_size)

    y[a1==1&s==1] <- stats::rbinom(length(which(a1==1&s==1)), size = 1, prob = response_prob[1])
    y[a1==1&s==0&a2==1] <- stats::rbinom(length(which(a1==1&s==0&a2==1)), size = 1, prob = response_prob[2])
    y[a1==1&s==0&a2==-1] <- stats::rbinom(length(which(a1==1&s==0&a2==-1)), size = 1, prob = response_prob[3])
    y[a1==-1&s==1] <- stats::rbinom(length(which(a1==-1&s==1)), size = 1, prob = response_prob[4])
    y[a1==-1&s==0&a2==1] <- stats::rbinom(length(which(a1==-1&s==0&a2==1)), size = 1, prob = response_prob[5])
    y[a1==-1&s==0&a2==-1] <- stats::rbinom(length(which(a1==-1&s==0&a2==-1)), size = 1, prob = response_prob[6])


  return(data.frame(a1,s,a2,y))
}

