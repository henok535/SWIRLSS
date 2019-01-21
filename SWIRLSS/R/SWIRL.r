
#' Sample Size Using SWIRL Method
#'
#'
#' This function is used to estimate the sample size needed to guarantee the
#' 95\% CI of the parameter theta is within used defined width.
#'
#'
#' @param modpar a vector containing the model coefficients
#' @param dist is for the distribution of the biomarker. Specify normal or uniform
#' @param mcdat is a vector of data with sample size arranged from the smallest to the largest
#' @param m number of data sets to generate for a given sample size
#' @param bmrkpar parameters for the distribution of the biomarker.
#' @param new is a vector for user specified width.
#' @return a numeric value of estimated sample size for a given used specificed width.
#' @export



ssswirl <- function(modpar, dist, mcdat, m, bmrkpar, new, ...) {
    
    B = list()
    finalresult = list()
    output = vector()
    
    # mcdat = seq(300,1500,100) # sample size starting from 300 to 1500 by an increament of 100
    
    
    for (i in 1:length(mcdat)) {
        
        set.seed(535)
        B[[i]] = logdata(m, mcdat[i], modpar, dist, bmrkpar)
        
        finalresult[[i]] = lapply(B[[i]], function(a) bootsmethod(a, dist, bmrkpar))
        
        output[i] <- c(colMeans(as.data.frame(matrix(unlist(finalresult[[i]]), ncol = 11, byrow = T))))[5]
        
    }
    
    
    width <- output
    w = 1/(width)^2
    predSamplesize = round(predict(lm(mcdat ~ w), new, level = 0.95, interval = "confidence"), digits = 0)
    
    return(predSamplesize)
    
}



