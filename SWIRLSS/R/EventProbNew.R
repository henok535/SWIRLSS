#' Probability of an Even Under a Given Treatment
#'
#'
#' This function is used to calculate the probability of an even under a given
#' treatment. If we set trt=0, that is if the default treatment is treat none,
#' if we set trt=1, for a default treatment 'Treat all' and if we set trt=3 for
#' the optimal treatment. We assumed here the biomarker to have a normal distribution
#' with mean zero and standard deviation one.
#'
#'
#' @param modpar is a vector having the value for the coefficients of the model
#' @param dist is for the distribution of the biomarker. Specify either normal or uniform
#' @param trtopt is for the treatment option. Specify treat_all, treat_none or treat_Opt.
#' @param bmrkpar parameters of the biomarker assuming a specific distribution
#' @return  a numberic probability value between 0 and 1
#' @export


eventprob <- function(modpar, dist, trtopt, bmrkpar, ...) {
    
    myreturn1 <- 0
    myretrun2 <- 0
    myreturn <- 0
    
    
    ####################################################################################### 
    
    # the functions to be integrated assuming a normal distribution
    
    # (a) when subjects are assigned to no-treatment
    
    fx1unif <- function(x) {
        
        a1 <- exp(modpar[1] + modpar[2] * x)
        a2 <- 1 + a1
        a1/a2
    }
    
    # (b) When subjects are assigned to treatment
    
    fx2unif <- function(x) {
        
        b1 <- exp(modpar[1] + modpar[3] + (modpar[2] + modpar[4]) * x)
        b2 <- 1 + b1
        b1/b2
    }
    
    ####################################################################################### 
    
    
    # (a) when subjects are assigned to no-treatment
    
    fx1norm <- function(x) {
        
        a1 <- exp(modpar[1] + modpar[2] * x)
        a2 <- 1 + a1
        
        b1 <- 1/(sqrt(2 * pi * bmrkpar[2]))
        b2 <- -(x - bmrkpar[1])^2
        b3 <- 2 * bmrkpar[2]
        
        (a1/a2) * b1 * exp(b2/b3)
        
    }
    
    # (b) When subjects are assigned to treatment
    
    fx2norm <- function(x) {
        
        a1 <- exp(modpar[1] + modpar[3] + (modpar[2] + modpar[4]) * x)
        a2 <- 1 + a1
        
        b1 <- 1/(sqrt(2 * pi * bmrkpar[2]))
        b2 <- -(x - bmrkpar[1])^2
        b3 <- 2 * bmrkpar[2]
        
        (a1/a2) * b1 * exp(b2/b3)
        
    }
    
    ################################################################################## 
    
    # We used simpson method to perform the infinite itegrals Function used for simpson integration
    
    simpson_v2 <- function(fun, a1, a2, n = 1000) {
        
        # assume a1 < a2 and n is an even positive integer
        
        if (a1 == -Inf & a2 == Inf) {
            
            f <- function(t) (fun((1 - t)/t) + fun((t - 1)/t))/t^2
            s <- simpson_v2(f, 0, 1, n)
            
        } else if (a1 == -Inf & a2 != Inf) {
            
            f <- function(t) fun(a2 - (1 - t)/t)/t^2
            s <- simpson_v2(f, 0, 1, n)
            
        } else if (a1 != -Inf & a2 == Inf) {
            
            f <- function(t) fun(a1 + (1 - t)/t)/t^2
            s <- simpson_v2(f, 0, 1, n)
        } else {
            
            h <- (a2 - a1)/n
            x <- seq(a1, a2, by = h)
            y <- fun(x)
            y[is.nan(y)] = 0
            s <- y[1] + y[n + 1] + 2 * sum(y[seq(2, n, by = 2)]) + 4 * sum(y[seq(3, n - 1, by = 2)])
            s <- s * h/3
            
        }
        
        return(s)
    }
    ####################################################################################### 
    
    # when the interaction coefficient is zero,it means that the biomarker is not useful for treatment
    # selection, so stop the program.
    
    if (modpar[4] == 0) {
        stop("Biomarker is not useful for treatment selection")
    }
    
    # Now calculate the probability of an event for a given treatment scenario
    
    # A) When the biomarker has uniform distribution
    
    # (1) treatment option is treat none
    
    if (dist == "uniform" & trtopt == "treat_none") {
        
        myreturn <- integrate(fx1unif, bmrkpar[1], bmrkpar[2])
        myreturn <- myreturn$value
        
        # (2) treatment option is treat all
        
    } else if (dist == "uniform" & trtopt == "treat_all") {
        
        myreturn <- integrate(fx2unif, bmrkpar[1], bmrkpar[2])
        myreturn <- myreturn$value
        
        # (3) treatment option is the optimal treatment rule this rule is set according the equations developed to
        # determine the optimal treatment assignment which also depends on the sign of the interaction coefficient
        # from the fitted logistic regression model.
        
    } else if (dist == "uniform" & trtopt == "treat_opt") {
        
        if (modpar[4] < 0) {
            
            c0 <- bmrkpar[1]
            d0 <- -modpar[3]/modpar[4]
            c1 <- -modpar[3]/modpar[4]
            d1 <- bmrkpar[2]
            
        } else {
            c0 <- -modpar[3]/modpar[4]
            d0 <- bmrkpar[2]
            c1 <- bmrkpar[1]
            d1 <- -modpar[3]/modpar[4]
        }
        
        
        myreturn1 <- integrate(fx2unif, lower = c1, upper = d1)
        myreturn2 <- integrate(fx1unif, lower = c0, upper = d0)
        myreturn <- myreturn1$value + myreturn2$value
        
        
        
        # B) When the biomarker has normal distribution
        
        # (4) treatment option is treat all
        
    } else if (dist == "normal" & trtopt == "treat_none") {
        
        myreturn <- simpson_v2(fx1norm, -Inf, Inf, n = 1000)
        myreturn <- myreturn
        
        
        
        # (5) treatment option is treat all
        
    } else if (dist == "normal" & trtopt == "treat_all") {
        
        myreturn <- simpson_v2(fx2norm, -Inf, Inf, n = 1000)
        myreturn <- myreturn
        
        # (6) treatment option is the optimal treatment rule this rule is set according the equations developed to
        # determine the optimal treatment assignment which also depends on the sign of the interaction coefficient
        # from the fitted logistic regression model.
        
    } else if (dist == "normal" & trtopt == "treat_opt") {
        
        if (modpar[4] < 0) {
            
            c0 <- -Inf
            d0 <- -modpar[3]/modpar[4]
            c1 <- -modpar[3]/modpar[4]
            d1 <- Inf
            
        } else {
            c0 <- -modpar[3]/modpar[4]
            d0 <- Inf
            c1 <- -Inf
            d1 <- -modpar[3]/modpar[4]
        }
        
        
        myreturn1 <- simpson_v2(fx2norm, c1, d1, 1000)
        myreturn2 <- simpson_v2(fx1norm, c0, d0, 1000)
        myreturn <- myreturn1 + myreturn2
        
    }
    
    
    return(myreturn)
}

############### End of the Code ###########################################################
