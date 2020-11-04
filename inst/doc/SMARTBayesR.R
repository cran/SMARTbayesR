## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(SMARTbayesR)

## -----------------------------------------------------------------------------
set.seed(23856)
dat <- SimDesign1(sample_size = 250,
           response_prob = c(0.2,0.3,0.4,0.5,0.6,0.7),
           stage_one_trt_one_response_prob = 0.7,
           stage_one_trt_two_response_prob = 0.4)

## -----------------------------------------------------------------------------
set.seed(39864)

posterior_trt_seq_draws <- PosteriorTrtSeqProb(niter = 10000,
                    dat,
                    design = "design-1")


## -----------------------------------------------------------------------------
posterior_EDTR_draws <- PosteriorEDTRProbs(posterior_trt_seq_draws)

## -----------------------------------------------------------------------------
MCBUpperLimits(thetadraws=posterior_EDTR_draws,
               alpha=0.05,
               design="design-1")

## -----------------------------------------------------------------------------
LogOR(response_prob = c(0.3,0.3,0.3,0.8,0.6,0.7),
      stage_one_trt_one_response_prob = 0.6,
      stage_one_trt_two_response_prob = 0.4,
      design = "design-1")

## -----------------------------------------------------------------------------
set.seed(2364)
PowerBayesian("design-1",
              sample_size = 100,
              response_prob = c(0.2,0.3,0.4,0.5,0.6,0.7),
              stage_one_trt_one_response_prob = 0.6,
              stage_one_trt_two_response_prob = 0.4,
              rejection_indices = c(1,2))

