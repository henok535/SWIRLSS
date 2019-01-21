#' Estimating the modparicients For the Logistic Regression MOdel
#'
#'
#' This Function is used to Generate the beta modparicients from the
#' user specified 'K' values. These Beta are then used to generate data
#' from the logistic regression model.
#'
#'
#' @param k  a vector of four user defined elements
#' @param z a vector of two elements, the 25th and 75th percetile value of the biomarker
#' @return  a vector of beta modparicients
#' @author Henok Woldu
#' @details
#' This function take the values of K1-k4 and Z1-Z2 and return the four beta values.
#' The Beta values correspond to the modparicients in the multiple logistic regression
#' we fitted. The first coef is for the intercept, the 2nd for the modpar of the biomarker
#' , the 3rd for the treatment and the last for the interaction term.
#' @export


coeffgen <- function(k, z) {
    
    Bvalues = numeric(0)
    
    
    if (k[1] == k[3] && k[2] == k[4]) {
        warning("Biomarker is not related to Outcome")
    }
    if (z[1] == z[2]) {
        stop("Wrong  choice of z values")
    }
    if (length(k)) 
        
    
    Bvalues[1] <- (k[2] * z[1] - k[1] * z[2])/(z[1] - z[2])
    Bvalues[2] <- (k[1] - k[2])/(z[1] - z[2])
    Bvalues[3] <- (k[1] * z[2] - k[2] * z[1] - k[3] * z[2] + k[4] * z[1])/(z[1] - z[2])
    Bvalues[4] <- (k[2] - k[1] + k[3] - k[4])/(z[1] - z[2])
    
    
    return(c(Bvalues))
}





#' Generating Data From Multiple Logistic Regression Model
#'
#'
#' This function is used to generate data from a multiple logistic
#' regression model using the modparicients obtained from the \code{Betas}
#' function.
#'
#'
#' @param m is the number of data set needed to be generated
#' @param n is the number of sample size needed
#' @param modpar is the vector of the beta model parameters (coefficients)
#' @param dist is the distribution of the biomarker
#' @param bmrkpar parameters of the biomarker distribution
#' @param bmrkdgn is the type of biomarker study design

#' @return  a list with 'm' data set each having 'n' sample size
#' @export


logdata <- function(m, n, modpar, dist, bmrkpar, bmrkdgn = "stratified", ...) {
    
    # dist : is for the distribution of the parameter. We have only considered here normal , gamma and uniform
    # distributions
    set.seed(535)
    X <- matrix(nrow = n, ncol = m)
    z <- matrix(nrow = n, ncol = m)
    outcome <- matrix(nrow = n, ncol = m)
    mylog <- matrix(nrow = n, ncol = m)
    A <- matrix(nrow = n, ncol = m)
    myreturn <- list()
    
    for (i in 1:m) {
        
        if (bmrkdgn == "strategy") {
            
            A[, i] <- sample(0:1, size = n, replace = T, prob = c(0.25, 0.75))
            
        } else if (bmrkdgn == "stratified") {
            
            A[, i] <- sample(rep(0:1, each = (n/2)), size = n, replace = F)
            
        }
        
        if (dist == "normal") {
            
            X[, i] <- rnorm(n, bmrkpar[1], bmrkpar[2])
            
        } else if (dist == "gamma") {
            
            X[, i] <- rgamma(n, bmrkpar[1], bmrkpar[2])
            
        } else if (dist == "uniform") {
            
            X[, i] <- runif(n, bmrkpar[1], bmrkpar[2])
        }
        
        z[, i] <- A[, i] * X[, i]  # for the biomarker by treatment interaction
        
        # calculate the probability of success p
        
        mylog[, i] <- ((exp(modpar[1] + modpar[2] * X[, i] + modpar[3] * A[, i] + modpar[4] * z[, i]))/(1 + 
            exp(modpar[1] + modpar[2] * X[, i] + modpar[3] * A[, i] + modpar[4] * z[, i])))
        
        # generate an oucome from a binomial distribution with success probability p
        outcome[, i] <- rbinom(n, 1, mylog[, i])
        
        myreturn[[i]] <- data.frame(outcome[, i], A[, i], X[, i], z[, i])
        
    }
    myreturn <- lapply(myreturn, setNames, nm = c("outcome", "A", "X", "Z"))
    
    
    return(myreturn)
    
}  # end of code


















