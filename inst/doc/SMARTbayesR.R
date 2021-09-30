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
               design="design-1",
               type="log-OR")

## -----------------------------------------------------------------------------
LogOR(response_prob = c(0.3,0.3,0.3,0.8,0.6,0.7),
      stage_one_trt_one_response_prob = 0.6,
      stage_one_trt_two_response_prob = 0.4,
      design = "design-1")

## ---- results="hide"----------------------------------------------------------
set.seed(2364)
power1 <- PowerBayesian("design-1",
              sample_size = 100,
              response_prob = c(0.2,0.3,0.4,0.5,0.6,0.7),
              stage_one_trt_one_response_prob = 0.6,
              stage_one_trt_two_response_prob = 0.4,
              type="log-OR",
              threshold=1)

## -----------------------------------------------------------------------------
power1

## -----------------------------------------------------------------------------
LogRR(response_prob = c(0.3,0.3,0.3,0.8,0.6,0.7),
      stage_one_trt_one_response_prob = 0.6,
      stage_one_trt_two_response_prob = 0.4,
      design = "design-1")

## ---- results="hide"----------------------------------------------------------
set.seed(23641)
power2 <- PowerBayesian("design-1",
              sample_size = 100,
              response_prob = c(0.2,0.3,0.4,0.5,0.6,0.7),
              stage_one_trt_one_response_prob = 0.6,
              stage_one_trt_two_response_prob = 0.4,
              type="log-RR",
              threshold=0.8)

## -----------------------------------------------------------------------------
power2

## -----------------------------------------------------------------------------
RD(response_prob = c(0.3,0.3,0.3,0.8,0.6,0.7),
      stage_one_trt_one_response_prob = 0.6,
      stage_one_trt_two_response_prob = 0.4,
      design = "design-1")

## ---- results="hide"----------------------------------------------------------
set.seed(236412)
power3 <- PowerBayesian("design-1",
              sample_size = 100,
              response_prob = c(0.2,0.3,0.4,0.5,0.6,0.7),
              stage_one_trt_one_response_prob = 0.6,
              stage_one_trt_two_response_prob = 0.4,
              type="RD",
              threshold=0.3)

## -----------------------------------------------------------------------------
power3

## -----------------------------------------------------------------------------
LogOR(response_prob = c(0.2,0.3,0.4,0.5,0.6,0.7,0.7,0.6),
      stage_one_trt_one_response_prob = 0.6,
      stage_one_trt_two_response_prob = 0.4,
      design = "general")

## ---- results="hide"----------------------------------------------------------
set.seed(23644)
power4 <- PowerBayesian("general",
              sample_size = 250,
              response_prob = c(0.2,0.3,0.4,0.5,0.6,0.7,0.7,0.6),
              stage_one_trt_one_response_prob = 0.6,
              stage_one_trt_two_response_prob = 0.4,
              type="log-OR",
              threshold=1)

## -----------------------------------------------------------------------------
power4

## -----------------------------------------------------------------------------
LogOR(response_prob = c(0.2,0.3,0.4,0.5,0.6,0.7,0.7,0.6,0.9),
      stage_one_trt_one_response_prob = 0.6,
      stage_one_trt_two_response_prob = 0.4,
      stage_one_trt_three_response_prob = 0.7,
      design = "design-3")

## ---- results="hide"----------------------------------------------------------
set.seed(236445)
power5 <- PowerBayesian("design-3",
              sample_size = 250,
              response_prob = c(0.2,0.3,0.4,0.5,0.6,0.7,0.7,0.6,0.9),
              stage_one_trt_one_response_prob = 0.6,
              stage_one_trt_two_response_prob = 0.4,
              stage_one_trt_three_response_prob = 0.7,
              type="log-OR",
              threshold=0.8)

## -----------------------------------------------------------------------------
power5

