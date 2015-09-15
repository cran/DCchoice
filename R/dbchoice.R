
dbchoice <- function(formula, data, dist = "log-logistic", par = NULL, ...){
  # evaluating the formula and stops if the formula is not defined correctly
  if (!inherits(formula, "formula")) stop("invalid formula")
  # stop if the LHS does not contain two variables
  if(length(formula[[2]]) != 3) stop("LHS variable in the formula must be like y1 + y2 ")
  
  # checking the distribution
  if(dist != "logistic" & dist != "log-logistic" & dist != "normal" & dist != "log-normal" & dist != "weibull"){
    stop("'dist' is incorrect.")
  }

  # extracting explanatory variables (including the intercept) from the specified data frame
  cl <- match.call()            # a call to the function
  if(missing(data)) stop("the name of the data frame object must be supplied in the 'data' argument")
  
  data <- eval(data, parent.frame())  # retrieving the data
  
  # removing observations with missing values
  na.num <- max(sum(as.numeric(is.na(data))))
  if(na.num != 0){ 
    d1 <- nrow(data)
    data <- na.omit(data)
    d2 <- nrow(data)
    warning(paste("Missing values detected.", d1 - d2, "rows are removed.", sep = " "))
  }

  # defining the dependent variable
  lhs1 <- formula[[2]][[2]]       # extracting from the formula the name of the variable for the yes/no to the first bid
  lhs2 <- formula[[2]][[3]]       # extracting from the formula the name of the variable for the yes/no to the second bid
  y1 <- eval(lhs1, data)          # yes/no to the first bid
  y2 <- eval(lhs2, data)          # yes/no to the second bid
  
  nobs <- length(y1)
  
  LP1 <- formula[[3]][[3]][[2]]   # extracting from the formula the name of the variable for the first bid
  LP2 <- formula[[3]][[3]][[3]]   # extracting from the formula the name of the variable for the second bid
  Lfirst <- eval(LP1, data)       # the first bids
  Lsecond <- eval(LP2, data)      # the second bids
  
  # making dummy variables for the first and second bids
  if(is.factor(y1)){   # when the yes/no variables are defined as factor
    yy <- ifelse(y1 == "yes" & y2 == "yes", 1, 0)
    yn <- ifelse(y1 == "yes" & y2 == "no", 1, 0)
    ny <- ifelse(y1 == "no" & y2 == "yes", 1, 0)
    nn <- ifelse(y1 == "no" & y2 == "no", 1, 0)
  } else {
    yy <- ifelse(y1 == 1 & y2 == 1, 1, 0)
    yn <- ifelse(y1 == 1 & y2 == 0, 1, 0)
    ny <- ifelse(y1 == 0 & y2 == 1, 1, 0)
    nn <- ifelse(y1 == 0 & y2 == 0, 1, 0)
   }

  BID <- ifelse(Lfirst > Lsecond, Lsecond, Lfirst)
    # }

    # Creating a design matrix
    bid <- cbind(Lfirst, Lsecond)   # the first and the second stage bids
    yvar <- cbind(yy, yn, ny, nn)   # yes/no to "bid"
    
    ff2 <- ~.
    ff2[[2]] <- formula[[3]][[2]]   # taking out the covariates from the formula
    if(ff2[[2]] != 1){ 
      X <- model.matrix(ff2, data)    # making a design matrix
    } else {
      X <- matrix(1, nrow = nobs, ncol = 1)
      colnames(X) <- "(Intercept)"
    }

  tmp.data <- data.frame(y1, X, BID)

   # obtaining initial parameter values by logit model
   if(is.null(par)){
        f.stage <- glm(y1~. -1, family = binomial(link = "probit"), data = tmp.data)
         ini <- f.stage$coefficients # saving initial values for ML estimation
         npar <- length(ini)
         ini[npar] <- ifelse(ini[npar] > 0, -ini[npar], ini[npar])     # gives a negative initial value for the bid coefficient
         if(substr(dist, 1, 4) == "log-") names(ini)[npar] <- "log(bid)"
         names(ini)[1] <- "(Intercept)"
   } else { # initial parameter values are supplied by the user
      if(length(par) != ncol(tmp.data)-1) stop("the length of 'par' does not coincide with the number of explanatory variables.")
      ini <- par
      f.stage <- ini
    }


  if(dist == "logistic" | dist == "log-logistic"){
      # likelihood function
      dbLL <- function(param, dvar, ivar, bid){
        yy <- dvar[, 1]
        yn <- dvar[, 2]
        ny <- dvar[, 3]
        nn <- dvar[, 4]

        X1 <- cbind(ivar, bid[, 1])
        X2 <- cbind(ivar, bid[, 2])
        
        ll <- 
        sum(plogis(-X2[yy == 1, ]%*%param, lower.tail = FALSE, log.p = TRUE))  + 
        sum(plogis(-X2[nn == 1, ]%*%param, lower.tail = TRUE, log.p = TRUE))   + 
        sum(log(plogis(-X2[yn == 1, ]%*%param, lower.tail = TRUE, log.p = FALSE) - 
                 plogis(-X1[yn == 1, ]%*%param, lower.tail = TRUE, log.p = FALSE))) +   
        sum(log(plogis(-X1[ny == 1, ]%*%param, lower.tail = TRUE, log.p = FALSE) - 
                 plogis(-X2[ny == 1, ]%*%param, lower.tail = TRUE, log.p = FALSE)))  
        ifelse(is.finite(ll), return(-ll), NaN) 
      }
  } else if(dist == "normal" | dist == "log-normal") {
      # likelihood function
        dbLL <- function(param, dvar, ivar, bid){
          yy <- dvar[, 1]
          yn <- dvar[, 2]
          ny <- dvar[, 3]
          nn <- dvar[, 4]

          X1 <- cbind(ivar, bid[, 1])
          X2 <- cbind(ivar, bid[, 2])
          
          ll <- 
          sum(pnorm(-X2[yy == 1, ]%*%param, lower.tail = FALSE, log.p = TRUE))  + 
          sum(pnorm(-X2[nn == 1, ]%*%param, lower.tail = TRUE, log.p = TRUE))   + 
          sum(log(pnorm(-X2[yn == 1, ]%*%param, lower.tail = TRUE, log.p = FALSE) - 
                   pnorm(-X1[yn == 1, ]%*%param, lower.tail = TRUE, log.p = FALSE))) +   
          sum(log(pnorm(-X1[ny == 1, ]%*%param, lower.tail = TRUE, log.p = FALSE) - 
                   pnorm(-X2[ny == 1, ]%*%param, lower.tail = TRUE, log.p = FALSE)))  
        ifelse(is.finite(ll), return(-ll), NaN) 
        }
  } else if(dist == "weibull"){
        dbLL <- function(param, dvar, ivar, bid){
          yy <- dvar[, 1]
          yn <- dvar[, 2]
          ny <- dvar[, 3]
          nn <- dvar[, 4]
          
          X1 <- cbind(ivar, bid[, 1])
          X2 <- cbind(ivar, bid[, 2])
          
          ll <- 
          sum(pweibull(exp(-X2[yy == 1, , drop = FALSE]%*%param), shape = 1, lower.tail = FALSE, log.p = TRUE))  + 
          sum(pweibull(exp(-X2[nn == 1, , drop = FALSE]%*%param), shape = 1, lower.tail = TRUE, log.p = TRUE))   + 
          sum(log(pweibull(exp(-X2[yn == 1, , drop = FALSE]%*%param), shape = 1, lower.tail = TRUE, log.p = FALSE) - 
                   pweibull(exp(-X1[yn == 1, , drop = FALSE]%*%param), shape = 1, lower.tail = TRUE, log.p = FALSE))) +   
          sum(log(pweibull(exp(-X1[ny == 1, , drop = FALSE]%*%param), shape = 1, lower.tail = TRUE, log.p = FALSE) - 
                   pweibull(exp(-X2[ny == 1, , drop = FALSE]%*%param), shape = 1, lower.tail = TRUE, log.p = FALSE)))  
        ifelse(is.finite(ll), return(-ll), NaN) 
        }
  }
  
  # ML estimation of double-bounded dichotomous choice
  suppressWarnings(
        dbchoice <- optim(ini, fn = dbLL, method="BFGS", hessian = TRUE, dvar = yvar, ivar = X, bid = bid, control=list(abstol=10^(-30)))
    )
  npar <- length(dbchoice$par)
  
  # arranging outcomes into a single list variable
   output <- list(
      f.stage = f.stage,            # the outcome of the initial estimation
      dbchoice  = dbchoice,         # the outcome from the optimization
      coefficients = dbchoice$par,  # the coefficient estimates
      call = cl,                    # the function call
      formula = formula,            # the defined model formula
      Hessian = dbchoice$hessian,   # the numerical Hessian at the estimates
      distribution = dist,          # the specified error distribution
      loglik = -dbchoice$value,     # log-likelihood at the estimates
      convergence = ifelse(dbchoice$convergence == 0, TRUE, FALSE),   # convergence status
      niter = dbchoice$counts,      # the number of iterations
      nobs = nobs,                  # the number of observations
      covariates = X,               # a matrix of the covariates
      bid = bid,                    # the suggested bid
      yn = cbind(y1, y2),           # the acceptance/rejection variable
      data.name = data              # the data matrix
      )

  class(output) <- "dbchoice"       # setting the object class
  return(output)
}

# a function for summarizing the outputs
summary.dbchoice <- function(object, ...){
  # obtaining necessary components from the object
  coef <- object$coefficients
  npar <- length(coef)
  se <- sqrt(diag(solve(object$Hessian)))   # standard errors of the estimates
  # object$nobs <- nrow(object$covariates)
  bid <- object$bid
  X <- object$covariates
  dist = object$distribution

  # estimating the null model
    formula_null <- object$formula
    formula_null[[3]][[2]] <- 1
    db_null <- dbchoice(formula_null, data = eval(object$data.name), dist = dist, par = coef[c(1, npar)])
  
  # function for obrtaining AIC and BIC
  akaike <- function(loglik, npar, k ){
    -2*loglik + k*npar
  }

  # computing various mean estimates for different error distributions
  if(dist == "log-logistic"){
    b <- coef[npar]         # the estimate on log(bid)
    Xb <- sum(colMeans(X)*coef[-npar])  # a matrix containing the estimates times the covariates
    func <- function(x) plogis(-(Xb + b*log(x)), lower.tail = FALSE)  # a function for integrals
    object$medianWTP <- exp(-Xb/b)
    object$meanWTP <- ifelse(abs(b) > 1, integrate(func, 0, Inf, stop.on.error = FALSE)$value, Inf)  # a numerical integral for finding the mean WTP
    object$trunc.meanWTP <- integrate(func, 0, exp(max(bid)))$value                                  # a numerical integral for finding the truncated mean WTP
    object$adj.trunc.meanWTP <- integrate(func, 0, exp(max(bid)))$value/plogis(-(Xb + b*max(bid)))   # a numerical integral for finding the truncated mean WTP with adjustment
  }  else if(dist == "log-normal"){
    b <- coef[npar]         # the estimate on log(bid)
    Xb <- sum(colMeans(X)*coef[-npar])
    func <- function(x) pnorm(-(Xb + b*log(x)), lower.tail = FALSE)
    object$medianWTP <- exp(-Xb/b)
    object$meanWTP <- integrate(func, 0, Inf, stop.on.error = FALSE)$value
    object$trunc.meanWTP <- integrate(func, 0, exp(max(bid)))$value
    object$adj.trunc.meanWTP <- integrate(func, 0, exp(max(bid)))$value/pnorm(-(Xb + b*max(bid)))
  } else if(dist == "logistic"){
    b <- coef[npar]         # the estimate on log(bid)
    Xb <- sum(colMeans(X)*coef[-npar])
    func <- function(x) plogis(-(Xb + b*x), lower.tail = FALSE)
    object$medianWTP <- -Xb/b
    object$meanWTP <- integrate(func, 0, Inf, stop.on.error = FALSE)$value
    object$trunc.meanWTP <- integrate(func, 0, max(bid))$value
    object$adj.trunc.meanWTP <- integrate(func, 0, max(bid))$value/plogis(-(Xb + b*max(bid)))
  }  else if(dist == "normal"){
    b <- coef[npar]         # the estimate on log(bid)
    Xb <- sum(colMeans(X)*coef[-npar])
    func <- function(x) pnorm(-(Xb + b*x), lower.tail = FALSE)
    object$medianWTP <- -Xb/b
    object$meanWTP <- integrate(func, 0, Inf)$value
    object$trunc.meanWTP <- integrate(func, 0, max(bid))$value
    object$adj.trunc.meanWTP <- integrate(func, 0, max(bid))$value/pnorm(-(Xb + b*max(bid)))
  } else if(dist == "weibull"){
    b <- coef[npar]
    Xb <- sum(colMeans(X[, -npar])*coef[-npar])
    
    func <- function(x) pweibull(exp(-Xb - b*log(x)), shape=1, lower.tail=FALSE)
    object$medianWTP <- exp(-Xb/b)*(log(2))^(-1/b)
    object$meanWTP <- ifelse(abs(b) > 1, exp(-Xb/b)*gamma(1-1/b), Inf)  
    object$trunc.meanWTP <- integrate(func, 0, exp(max(bid)), stop.on.error = FALSE)$value                         
    object$adj.trunc.meanWTP <- integrate(func, 0, exp(max(bid)), stop.on.error = FALSE)$value/pweibull(exp(-Xb - b*max(bid)), shape=1)                    
  }

  # computing pseudo-R^2
#   object$psdR2 <- 1 - object$loglik/db_null$loglik
#   names(object$psdR2) <- "pseudo-R^2 measure"
#   object$adjpsdR2 <- 1 - (object$loglik - npar)/db_null$loglik
#   names(object$adjpsdR2) <- "adjusted pseudo-R^2 measure"
  
  # Likelihood Ratio Statistic
  LR <- -2*(db_null$loglik - object$loglik)
  d.f <- length(object$coefficients) - length(db_null$coefficients)
  pvalLR <- pchisq(LR, df = d.f, lower.tail = FALSE)
  object$LR.test <- c(LR, d.f, pvalLR)
  
  # creating a table for coefficients, se, etc. 
  zstat <- coef/se  # z statistics
  pval <- round(2*pnorm(-abs(zstat)), 6)  # p-value
  coef <- cbind(coef, se, zstat, pval)
  colnames(coef) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  
  object$coef <- coef
  
  # computing AIC and BIC by the function AKAIKE
  object$AIC <- akaike(object$loglik, npar, k = c(2, log(object$nobs)))
  names(object$AIC) <- c("AIC", "BIC")  
  
  class(object) <- "summary.dbchoice"
  return(object)
  
}

print.dbchoice <- function(x, digits = max(3, getOption("digits") - 1), ...){
  if(!x$convergence)  cat("The optimization did not converge\n")
  cat("\nDistribution:", x$distribution, "\n", sep = " ")
  print.default(format(x$coef, digits = digits), print.gap = 1, quote = FALSE)

  invisible(x)
}

print.summary.dbchoice <- function(x, digits = max(3, getOption("digits") - 1), ...){

  # what shall we do for Pearson residuals?

  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  cat("Formula:", deparse(x$formula, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")

  if(!x$convergence)  cat("The optimization did not converge\n")
  cat("Coefficients:", "\n")
  #print.default(format(x$coef, digits = digits), print.gap = 2, quote = FALSE, right = TRUE)
  printCoefmat(x$coef, digits = 4, dig.tst = 4)    # printCoefmat() is defined in stats package
  
  cat("\nDistribution:", x$distribution, "", sep = " ")
  cat("\nNumber of Obs.:", formatC(x$nobs, digits = 0), "\n")
  cat("Log-likelihood:", formatC(x$loglik, format="f", digits = digits), "\n")
#   cat("pseudo-R^2:", formatC(x$psdR2, format="f", digits = 4), 
#       ", adjusted pseudo-R^2:", formatC(x$adjpsdR2, format="f", digits = 4), "")
  cat("\nLR statistic:", formatC(x$LR.test[1], format="f", digits = 3), "on", formatC(x$LR.test[2], digits = 0), 
    "DF, p-value:", formatC(x$LR.test[3], format="f", digits = 3), "\n")
  cat("AIC:", formatC(x$AIC[1], format="f", digits = digits), ", BIC:", formatC(x$AIC[2], format="f", digits = digits), "\n")
  cat("\nIterations:", formatC(x$niter, digits = 0), "")
  cat("\nConvergence:", x$convergence, "\n")

  cat("\nWTP estimates:")
  if(is.finite(x$meanWTP)){
    cat("\n Mean :", x$meanWTP, "", sep = " ")
  } else {
    cat("\n Mean :", x$meanWTP, "(because of |beta_Lbid| < 1)", sep = " ")
  }
  cat("\n Mean :", x$trunc.meanWTP, "(truncated at the maximum bid)", sep = " ")
  cat("\n Mean :", x$adj.trunc.meanWTP, "(truncated at the maximum bid with adjustment)", sep = " ")
  cat("\n Median:", x$medianWTP, "\n", sep = " ")

}

######################################################################
