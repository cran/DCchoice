# a binary choice model for Single-bounded data. a simple logit or probit model
sbchoice <- function(formula, data, dist = "log-logistic", ...){
  # evaluating the formula and stops if the formula is not defined correctly
  if (!inherits(formula, "formula")) stop("invalid formula")

  cl <- match.call()            
#  if(missing(data)) data <- environment(formula)
  if(missing(data)) stop("the name of the data frame object must be supplied in the 'data' argument")
  
  data <- eval(data, parent.frame())

  # removing observations with missing values
  na.num <- max(sum(as.numeric(is.na(data))))
  if(na.num != 0){ 
    d1 <- nrow(data)
    data <- na.omit(data)
    d2 <- nrow(data)
    warning(paste("Missing values detected.", d1 - d2, "rows are removed.", sep = " "))
  }

  # defining the dependent variable 
  lhs <- formula[[2]]       # extracting the name of the acceptance/rejection variable from the formula supplied
  y1 <- eval(lhs, data)     # yes/no to the bids

  nobs <- length(y1)
  
  LP1 <- formula[[3]][[3]]  # extracting the name of the bid variable from the formula supplied
  BID <- eval(LP1, data)    # the suggested bid
  
  # handling the data matrix
  f2 <- ~.
  f2[[2]] <- formula[[3]][[2]]   # retrieving the names of the covariates other than the bid variable
  if(f2[[2]] != 1){
    X <- model.matrix(f2, data)    # making a design matrix
  } else {                         # for a model without covariate
    X <- matrix(1, nrow = nobs, ncol = 1) # making a constant term
    colnames(X) <- "(Intercept)"  # naming the constant
  }

  tmp.data <- data.frame(y1, X, BID)  # re-combining all the ables into a single dataframe

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
if(dist == "logistic" | dist == "log-logistic"){  # logistic or log-logistic error distribution
    glm.out <- glm(y1~. -1, family = binomial(link = "logit"), data = tmp.data)  # unrestricted model
    glm.null <- update(glm.out, .~ 1)  # restricted (only the intercept) model
    
    npar <- length(glm.out$coefficients) # the number of parameters
    if(substr(dist, 1, 4) == "log-") names(glm.out$coefficients)[npar] <- "log(bid)" # changing the name label if the model is log-logistic
    names(glm.out$coefficients)[1] <- "(Intercept)"
    estimates <- glm.out$coefficients  # taking out the estimates
    
} else if(dist == "normal" | dist == "log-normal") {  # normal or log-normal error distribution
    glm.out <- glm(y1~. -1, family = binomial(link = "probit"), data = tmp.data)
    glm.null <- update(glm.out, .~ 1)
    
    npar <- length(glm.out$coefficients)
    if(substr(dist, 1, 4) == "log-") names(glm.out$coefficients)[npar] <- "log(bid)"
    names(glm.out$coefficients)[1] <- "(Intercept)"
    estimates <- glm.out$coefficients
} else if(dist == "weibull"){
    # likelihood function
    sbLL <- function(param, dvar, ivar){
          y1 <- dvar
          X <- ivar
          ll <- 
            sum(pweibull(exp(-X[y1==1, , drop=FALSE]%*%param), shape = 1, lower.tail = FALSE, log.p=TRUE)) + 
            sum(pweibull(exp(-X[y1==0, , drop=FALSE]%*%param), shape = 1, lower.tail = TRUE, log.p=TRUE))
          ifelse(is.finite(ll), return(-ll), NaN) 
      }
    # initial parameter values
    ini.par <- glm(y1~. -1, family = binomial(link = "probit"), data = tmp.data)$coefficients
    ini.par.null <- glm(y1 ~ 1, family = binomial(link = "probit"), data = tmp.data)$coefficients
    
    # ML estimation
    suppressWarnings( # "glm." is nothing to do with GLM. The naming is merely because of compatibility
#        glm.out <- optim(ini.par, fn = sbLL, method="BFGS", hessian = TRUE, dvar = y1, ivar = cbind(X, BID), control=list(abstol=10^(-30)))
        glm.out <- optim(ini.par, fn = sbLL, method="BFGS", hessian = TRUE, dvar = y1, ivar = cbind(X, BID))
    )
    suppressWarnings(
#        glm.null <- optim(ini.par.null, fn = sbLL, method = "BFGS", hessian = TRUE, dvar = y1, ivar = matrix(1, nobs, 1), control=list(abstol=10^(-30)))
        glm.null <- optim(ini.par.null, fn = sbLL, method = "BFGS", hessian = TRUE, dvar = y1, ivar = matrix(1, nobs, 1))
    )
    names(glm.out$par)[1] <- "(Intercept)"
    if(dist == "weibull") {
      npar <- length(glm.out$par)     # the number of parameters
      names(glm.out$par)[npar] <- "log(bid)"
    }
    estimates <- glm.out$par
    # compatibility with other distributions
    glm.out$converged <- ifelse(glm.out$convergence == 0, TRUE, FALSE)  # convergence status of the ML estimation
    glm.out$iter <- glm.out$counts
    glm.out$coefficients <- glm.out$par
    glm.null$coefficients <- glm.null$par
} else {
    stop("dist must be logistic, normal or weibull")
}

  # arranging outcomes into a single list variable
   output <- list(
      coefficients = estimates, # the coefficient estimates
      call = cl,            # the function call
      formula = formula,    # the defined model formula
      glm.out = glm.out,    # the outcome of the unrestricted model
      glm.null = glm.null,  # the outcome of the null model
      distribution = dist,  # the specified error distribution
      covariates = X,       # a matrix of the covariates
      bid = BID,            # the suggested bid
      nobs = nobs,          # the number of observations
      yn = y1,              # the acceptance/rejection variable
      data.name = data      # the data matrix
      )

  class(output) <- "sbchoice"   # setting the object class
  return(output)
}

# summarizing the output
summary.sbchoice <- function(object, ...){
  dist <- object$distribution
  coef <- object$coefficients
  npar <- length(coef)
  X <- object$covariates
  bid <- object$bid

    # function for obrtaining AIC and BIC
    akaike <- function(loglik, npar, k ){
      -2*loglik + k*npar
    }

  if(dist == "weibull"){
    # creating a table for coefficients, se, etc. 
    se <- sqrt(diag(solve(object$glm.out$hessian)))
    zstat <- coef/se  # z statistics
    pval <- round(2*pnorm(-abs(zstat)), 6)  # p-value
    coefmat <- cbind(coef, se, zstat, pval)
    colnames(coefmat) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    object$coefmat <- coefmat
    
    object$loglik <- -object$glm.out$value
    object$loglik.null <- -object$glm.null$value

    
    b <- coef[npar]
    Xb <- sum(colMeans(X)*coef[-npar])
    
    func <- function(x) pweibull(exp(-Xb - b*log(x)), shape=1, lower.tail=FALSE)
    object$medianWTP <- exp(-Xb/b)*(log(2))^(-1/b)
    object$meanWTP <- ifelse(abs(b) > 1, exp(-Xb/b)*gamma(1-1/b), Inf)  
    object$trunc.meanWTP <- integrate(func, 0, exp(max(bid)), stop.on.error = FALSE)$value                         
    object$adj.trunc.meanWTP <- integrate(func, 0, exp(max(bid)), stop.on.error = FALSE)$value/pweibull(exp(-Xb - b*max(bid)), shape=1)                    
  } else {
    # obtaining necessary components from the object
    object$glm.summary <- summary.glm(object$glm.out)
    object$glm.null.summary <- summary.glm(object$glm.null)
    object$loglik <- logLik(object$glm.out)[1]
    object$loglik.null <- logLik(object$glm.null)[1]

    # computing various mean estimates for different error distributions
    if(dist == "log-logistic"){
      b <- coef[npar]         # the estimate on log(bid)
      Xb <- sum(colMeans(X[, -npar, drop = FALSE])*coef[-npar]) # a matrix containing the estimates times the covariates
      func <- function(x) plogis(-(Xb + b*log(x)), lower.tail = FALSE)  # a function for integrals
      object$medianWTP <- exp(-Xb/b)
      object$meanWTP <- ifelse(abs(b) > 1, integrate(func, 0, Inf, stop.on.error = FALSE)$value, Inf)   # a numerical integral for finding the mean WTP
      object$trunc.meanWTP <- integrate(func, 0, exp(max(bid)), stop.on.error = FALSE)$value            # a numerical integral for finding the truncated mean WTP
      object$adj.trunc.meanWTP <- integrate(func, 0, exp(max(bid)), stop.on.error = FALSE)$value/plogis(-(Xb + b*max(bid))) # a numerical integral for finding the truncated mean WTP with adjustment
    }  else if(dist == "log-normal"){
      b <- coef[npar]         # the estimate on log(bid)
      Xb <- sum(colMeans(X[, -npar, drop = FALSE])*coef[-npar])
      func <- function(x) pnorm(-(Xb + b*log(x)), lower.tail = FALSE)
      object$medianWTP <- exp(-Xb/b)
      object$meanWTP <- integrate(func, 0, Inf, stop.on.error = FALSE)$value
      object$trunc.meanWTP <- integrate(func, 0, exp(max(bid)), stop.on.error = FALSE)$value
      object$adj.trunc.meanWTP <- integrate(func, 0, exp(max(bid)), stop.on.error = FALSE)$value/pnorm(-(Xb + b*max(bid)))
    } else if(dist == "logistic"){
      b <- coef[npar]         # the estimate on log(bid)
      Xb <- sum(colMeans(X[, -npar, drop = FALSE])*coef[-npar])
      func <- function(x) plogis(-(Xb + b*x), lower.tail = FALSE)
      object$medianWTP <- -Xb/b
      object$meanWTP <- integrate(func, 0, Inf, stop.on.error = FALSE)$value
      object$trunc.meanWTP <- integrate(func, 0, max(bid))$value
      object$adj.trunc.meanWTP <- integrate(func, 0, max(bid))$value/plogis(-(Xb + b*max(bid)))
    }  else if(dist == "normal"){
      b <- coef[npar]         # the estimate on log(bid)
      Xb <- sum(colMeans(X[, -npar, drop = FALSE])*coef[-npar])
      func <- function(x) pnorm(-(Xb + b*x), lower.tail = FALSE)
      object$medianWTP <- -Xb/b
      object$meanWTP <- integrate(func, 0, Inf, stop.on.error = FALSE)$value
      object$trunc.meanWTP <- integrate(func, 0, max(bid), stop.on.error = FALSE)$value
      object$adj.trunc.meanWTP <- integrate(func, 0, max(bid), stop.on.error = FALSE)$value/pnorm(-(Xb + b*max(bid)))
    }
  } 
    # computing pseudo-R^2
    object$psdR2 <- 1 - object$loglik/object$loglik.null
    names(object$psdR2) <- "pseudo-R^2 measure"
    object$adjpsdR2 <- 1 - (object$loglik - npar)/object$loglik.null
    names(object$adjpsdR2) <- "adjusted pseudo-R^2 measure"

    # computing Likelihood Ratio Statistic and its p-value
    LR <- 2*(object$loglik - object$loglik.null)    # the LR statistic
    d.f <- length(object$glm.out$coefficients) - length(object$glm.null$coefficients) # the degrees of freedom
    pvalLR <- pchisq(LR, df = d.f, lower.tail = FALSE)  # p-value
    object$LR.test <- c(LR, d.f, pvalLR)

    # computing AIC and BIC by the function AKAIKE
    object$AIC <- akaike(object$loglik, npar, k = c(2, log(object$nobs)))
    names(object$AIC) <- c("AIC", "BIC")  

    class(object) <- "summary.sbchoice"
    return(object)
}


print.sbchoice <- function(x, digits = max(3, getOption("digits") - 1), ...){
  cat("\nDistribution:", x$distribution, "\n", sep = " ")
#  print(x$glm.out$coefficients)
  print(x$coefficients)
  invisible(x)
}


print.summary.sbchoice <- function(x, digits = max(3, getOption("digits") - 1), ...){
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  cat("Formula:", deparse(x$formula, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  
  if(!x$glm.out$converged)  cat("The optimization did not converge\n")
  cat("Coefficients:", "\n")
  if(x$distribution == "weibull" | x$distribution == "log-weibull") {
    printCoefmat(x$coefmat, digits = 4, dig.tst = 3)
  } else {
    printCoefmat(x$glm.summary$coefficients, digits = 4, dig.tst = 3)    # printCoefmat() is defined in stats package
  }
  cat("\nDistribution:", x$distribution, "", sep = " ")
  cat("\nNumber of Obs.:", formatC(x$nobs, digits = 0), " ")
  cat("\nlog-likelihood:", x$loglik, "\n", sep = " ")
  cat("pseudo-R^2:", formatC(x$psdR2, format="f", digits = 4), 
      ", adjusted pseudo-R^2:", formatC(x$adjpsdR2, format="f", digits = 4), "")
  cat("\nLR statistic:", round(x$LR.test[1], 3), "on", formatC(x$LR.test[2], digits = 0), 
    "DF, p-value:", formatC(x$LR.test[3], format="f", digits = 3), "\n")
  cat("AIC:", formatC(x$AIC[1], format="f", digits = digits), ", BIC:", formatC(x$AIC[2], format="f", digits = digits), "\n")
  cat("\nIterations:", x$glm.out$iter, " ")
  cat("\nConvergence:", x$glm.out$converged, "\n")
  
  cat("\nWTP estimates:")
  if(is.finite(x$meanWTP)){
    cat("\n Mean :", x$meanWTP, "", sep = " ")
  } else {
    cat("\n Mean :", x$meanWTP, "(because of |beta_bid| < 1)", sep = " ")
  }
  cat("\n Mean :", x$trunc.meanWTP, "(truncated at the maximum bid)", sep = " ")
  cat("\n Mean :", x$adj.trunc.meanWTP, "(truncated at the maximum bid with adjustment)", sep = " ")
  cat("\n Median :", x$medianWTP, "\n", sep = " ")

}


