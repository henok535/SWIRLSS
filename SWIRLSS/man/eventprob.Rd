% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EventProbNew.R
\name{eventprob}
\alias{eventprob}
\title{Probability of an Even Under a Given Treatment}
\usage{
eventprob(modpar, dist, trtopt, bmrkpar, ...)
}
\arguments{
\item{modpar}{is a vector having the value for the coefficients of the model}

\item{dist}{is for the distribution of the biomarker. Specify either normal or uniform}

\item{trtopt}{is for the treatment option. Specify treat_all, treat_none or treat_Opt.}

\item{bmrkpar}{parameters of the biomarker assuming a specific distribution}
}
\value{
a numberic probability value between 0 and 1
}
\description{
This function is used to calculate the probability of an even under a given
treatment. If we set trt=0, that is if the default treatment is treat none,
if we set trt=1, for a default treatment 'Treat all' and if we set trt=3 for
the optimal treatment. We assumed here the biomarker to have a normal distribution
with mean zero and standard deviation one.
}
