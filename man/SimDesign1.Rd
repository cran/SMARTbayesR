% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SimDesign1.R
\name{SimDesign1}
\alias{SimDesign1}
\title{Simulate a SMART with Design 1}
\usage{
SimDesign1(
  sample_size = 250,
  response_prob = c(0.5, 0.5, 0.5, 0.8, 0.7, 0.5),
  stage_one_trt_one_response_prob = 0.7,
  stage_one_trt_two_response_prob = 0.4
)
}
\arguments{
\item{sample_size}{the total sample size}

\item{response_prob}{a vector of probabilities of response for each of the 6 embedded treatment sequences.}

\item{stage_one_trt_one_response_prob}{the probability of response to first stage-1 treatment.}

\item{stage_one_trt_two_response_prob}{the probability of response to second stage-1 treatment.}
}
\value{
a data frame with treatment response indicators for each stage of treatment, a1 and a2,
end of stage-1 response indicator s, and final outcome y.
}
\description{
This function simulates a SMART with Design 1.
}
\examples{
dat <- SimDesign1(sample_size=250,
                             response_prob = c(0.5,0.9,0.3,0.7,0.5,0.8),
                             stage_one_trt_one_response_prob = 0.7,
                             stage_one_trt_two_response_prob = 0.4)

}
