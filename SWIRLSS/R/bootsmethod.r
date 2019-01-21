
#' The SWIRL Sample Size Estimation Method
#'
#'
#' This function is used to estimate the sample size that guarantee the
#' 95\% CI of the paramter Theta to be within user defined width. It first
#' generate a data set from the model used to calculate Theta with a sample size ranging from
#' 300 to 3500 and calculate the 95|% CI width of each sample. Using OLS , then a model is fitted
#' with n as response and W as a covariate.
#'
#'
#' @param dist is for the distribution of the biomarker. Specify either normal or uniform
#' @param bmrkpar is a vector containing the parameters of the distribution.
#' @param data is a data frame containing the data set
#' @return  a list containing values for Theta, 95\% interval for Theta and the converage probability for Theta
#' @export


#SWIRL sample Size Method


bootsmethod <- function(data, dist, bmrkpar,...) {

       getcoef = function(data,indeces) {

       # fit a logistic regression model and extract
       # the coefficients from the fitted model

         dat = data[indeces,]
         fit = glm(dat$outcome~dat$X+dat$A+dat$Z, family = binomial(link = "logit"))
         betas = c(fit$coef);

    # for a given biomarker distribution estimate the probability of a unfavorable
    # outcome under each treatment condition : assuming the dafault is treat none, treal all and
    # under the optimal treatment ( biomarker guided treatment)

           if (dist == "normal") {

                       A0 = eventprob(betas, dist = "normal", trt = "treat_none", bmrkpar);
                       A1 = eventprob(betas, dist = "normal", trt = "treat_all", bmrkpar);
                      Aopt = eventprob(betas, dist = "normal", trt = "treat_opt", bmrkpar);

                      Theta1 = A1 - Aopt  # estimate Theta assuming the default treatment is treat all
                      Theta0 = A0 - Aopt  # estimate Theta assuming the default treatment is treat all
                      output = c(A0,A1,Aopt, Theta1, Theta0)

            } else if (dist == "uniform") {

                         A0 = eventprob(betas, dist = "uniform", trt = "treat_none", bmrkpar);
                         A1 = eventprob(betas, dist = "uniform", trt = "treat_all", bmrkpar);
                        Aopt = eventprob(betas, dist = "uniform", trt = "treat_opt", bmrkpar);

                  Theta1 = A1 - Aopt  # estimate Theta assuming the default treatment is treat all
                  Theta0 = A0 - Aopt  # estimate Theta assuming the default treatment is treat all
                  output = c(A0,A1,Aopt, Theta1, Theta0)

          }

         return(output)

  }

  # use boots trap method in order to get a 95% CI for Theta1 and Theta0

            set.seed(535)
          one.boot = boot(data, getcoef, R = 500)
          margProb = as.data.frame(one.boot$t)
            MeanA0 = mean(margProb$V1)
            MeanA1 = mean(margProb$V2)
          MeanAopt = mean(margProb$V3)
         MeanTheta1 = mean(margProb$V4)
         MeanTheta0 = mean(margProb$V5)


     coveragetheta1 = boot.ci(one.boot, conf = 0.95, type = "norm", index = 4)
     coveragetheta0 = boot.ci(one.boot, conf = 0.95, type = "norm", index = 5)


                Normalcovtheta1 = c(coveragetheta1[[4]][2],coveragetheta1[[4]][3])
                Normalcovtheta0 = c(coveragetheta0[[4]][2], coveragetheta0[[4]][3])
                    width1 = mean(Normalcovtheta1[2]) - mean(Normalcovtheta1[1])
                    width0 = mean(Normalcovtheta0[2]) - mean(Normalcovtheta0[1])

              finalresult = list( Treatnone = MeanA0, Treatall = MeanA1, OptimalTreat = MeanAopt,
                                 Theta1 = MeanTheta1, Widtheta1 = width1,
                                 CI4theta1 = c(mean(Normalcovtheta1[1]), mean(Normalcovtheta1[2])),
                                 Theta0 = MeanTheta0, Widtheta0 = width0,
                                 CI4theta0 = c(mean(Normalcovtheta0[1]), mean(Normalcovtheta0[2]))

                              )
     finalresult;

}   # End of code








