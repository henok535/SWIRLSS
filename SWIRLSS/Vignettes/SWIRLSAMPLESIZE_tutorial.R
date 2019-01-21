## ---- eval=FALSE , echo=T------------------------------------------------
#  
#  install.packages("SWIRLSS")  # install it first
#  
#  library(SWIRLSS)             # load the package
#  

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  
#  install.packages("devtools")  # install the devtools package if you have not done it before.
#  
#  library(devtools)            # load the package devtools
#  
#  # to install the SWIRLSS from my github acount ("henok535/SWIRLSS")
#  
#  install_github("henok535/SWIRLSS")
#  
#  library(SWIRLSS)
#  
#  # Well come to SWIRLSS package. It is now ready for use!
#  

## ------------------------------------------------------------------------

library(boot)


## ------------------------------------------------------------------------

# Here are some possible user defined values ( 11, k2, k3 and k3). These can be adjusted to meet specific
# problem under consideration.

kval1 <- c(K1 = log(0.25/0.75), k2 = log(0.75/0.25), k3 = log(0.75/0.25), k4= log(0.25/0.75))
kval2 <- c(k1 = log(0.10/0.90), k2 = log(0.90/0.10), k3 = log(0.90/0.10), k4 = log(0.10/0.90))
kval3 <- c(k1 = log(0.10/0.90), k2 = log(0.55/0.45), k3 = log(0.90/0.10), k4 = log(0.45/0.55))
kval4 <- c(k1 = log(0.25/0.75), k2 = log(0.75/0.25), k3 = log(0.50/0.50), k4 = log(0.50/0.50))
kval5 <- c(k1 = log(0.90/0.10), k2 = log(0.45/0.55), k3 = log(0.10/0.90), k4 = log(0.55/0.45))
kval6 <- c(k1 = log(0.60/0.40), k2 = log(0.40/0.60), k3 = log(0.50/0.50), k4 = log(0.50/0.50))


# get the 25th and 75th percentile value of the biomarker value  assuming it is
# is normaly or uniformly distributed

set.seed(535)
pardist <- c(0,1)    # standard uniform / normal. We can change this one as needed

bmrkunif <- runif(100000,pardist[1],pardist[2]) # Assuming the biomarker has uniform distribution
bmrknorm <- rnorm(100000,pardist[1],pardist[2]) # Assuming the biomarker has normal distribution

zval <- c(quantile(bmrkunif)[2], quantile(bmrkunif)[4]) #25th and 75th percentile biomarker

modcoeff <- coeffgen(kval1, zval)
modcoeff    # this provideds the coefficients beta0, beta1, beta2 and beta3
            # which will be used to generate data from the logistic model specified above.

## ------------------------------------------------------------------------

# When the biomarker distribution is uniform

Event_None <- eventprob(modpar = modcoeff, dist = "uniform", trtopt = "treat_none",  bmrkpar = pardist)
Event_All <- eventprob(modpar = modcoeff, dist = "uniform", trtopt = "treat_all",  bmrkpar = pardist)
Event_Opt <- eventprob(modpar = modcoeff, dist = "uniform",  trtopt = "treat_opt",  bmrkpar = pardist)

# then Theta1 and Theta0 will be

Theta1_true <- Event_All - Event_Opt
Theta0_true <- Event_None - Event_Opt
c(Theta1_true, Theta0_true)

# When the biomarker is generated from a normal distribution, one can uncomment the following codes
# to calculate the event probability under each scenario the their respective thetas.

  #Event_None <- eventprob(modpar = modcoeff, dist = "normal", trtopt = "treat_none",  bmrkpar = pardist)
  #Event_All <- eventprob(modpar = modcoeff, dist = "normal", trtopt = "treat_all",  bmrkpar = pardist)
  #Event_Opt <- eventprob(modpar = modcoeff, dist = "normal", trtopt = "treat_opt",  bmrkpar = pardist)

# then Theta1 and Theta0 will be

    #Theta1_true <- Event_All - Event_Opt
    #Theta0_true <- Event_None - Event_Opt
    #c(Theta1_true, Theta0_true)


## ------------------------------------------------------------------------

m = 25      # number of data sets
n = 500     # sample size for each data set


mydat <- logdata(m, n, modcoeff, "uniform", pardist)

datt <- mydat[[1]]     # sample of one data set            


## ------------------------------------------------------------------------

# alldatout <- lapply(mydat, function(a) bootsmethod(a, "uniform", pardist))
# outdat <- as.data.frame(matrix(unlist(alldatout), ncol = 11, byrow = T))
# 
# Finaloutput <- list(Treat_None = mean(outdat$V1),
#                     Treat_All = mean(outdat$V2),
#                     Treat_Opt = mean(outdat$V3),
#                     Theta1 = mean(outdat$V4),
#                     Theta195CI = c(mean(outdat$V6), mean(outdat$V7)),
#                     width_Theta1 = mean(outdat$V5),
#                     Converage_Theta1 = length(intersect(which(outdat$V6 < Theta1_true),
#                                                 which(outdat$V7 > Theta1_true))),
#                     Theta0 = mean(outdat$V8),
#                     Theta095CI = c(mean(outdat$V10), mean(outdat$V11)),
#                     width_Theta0 = mean(outdat$V9),
#                     Converage_Theta0 = length(intersect(which(outdat$V10 < Theta0_true),
#                                             which(outdat$V11 > Theta0_true)))
#                   )
# 
# Finaloutput



## ------------------------------------------------------------------------

# newdat1 <- data.frame(w = c(0.10,0.15,0.20)^-2)
# mcdat = seq(100,1500,50)
# samplesize1 <- ssswirl(modcoeff, "uniform", mcdat, m = 10, pardist,
#                                    bmrkdgn = "strategy", newdat1)
# 
# 
# samplesize1


