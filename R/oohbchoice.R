oohbchoice <- function (formula, data, subset, na.action = na.omit, dist = "log-logistic", par = NULL, ...){
# argument "na.action" was added in June 2016

  if (!inherits(formula, "Formula"))
    formula <- Formula(formula)

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
  
  mf <- match.call(expand.dots = TRUE)

# Revised in June 2016
#  m <- match(c("formula", "data", "subset"), names(mf), 0L)
m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)

  mf <- mf[c(1L, m)]
  mf$formula <- formula
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  original.data <- data
  data <- mf
  
#  # removing observations with missing values
#  na.num <- max(sum(as.numeric(is.na(data))))
#  if(na.num != 0){ 
#    d1 <- nrow(data)
#    data <- na.omit(data)
#    d2 <- nrow(data)
#    warning(paste("Missing values detected.", d1 - d2, "rows are removed.", sep = " "))
#  }

  # defining the dependent variable
  y1 <- model.part(formula, data = data, lhs = 1)[[1]]  # yes/no to the first bid
  y2 <- model.part(formula, data = data, lhs = 1)[[2]]  # yes/no to the second bid
  
  nobs <- length(y1)

Llow <- model.matrix(formula, data = data, rhs = 2)[, 2]
Lhig <- model.matrix(formula, data = data, rhs = 2)[, 3]
  
#  # making dummy variables for the first and second bids
#  if(is.factor(y1)){   # when the yes/no variables are defined as factor
#    yy <- ifelse(y1 == "yes" & y2 == "yes", 1, 0)
#    yn <- ifelse(y1 == "yes" & y2 == "no", 1, 0)
#    ny <- ifelse(y1 == "no" & y2 == "yes", 1, 0)
#    nn <- ifelse(y1 == "no" & y2 == "no", 1, 0)
#  } else {
#    yy <- ifelse(y1 == 1 & y2 == 1, 1, 0)
#    yn <- ifelse(y1 == 1 & y2 == 0, 1, 0)
#    ny <- ifelse(y1 == 0 & y2 == 1, 1, 0)
#    nn <- ifelse(y1 == 0 & y2 == 0, 1, 0)
#   }

if (is.factor(y1)) {
  # Respondents who face a lower bid in the first question
  firstLower.yy  <- ifelse(y1 == "yes" & y2 == "yes",  1, 0)
  firstLower.yn  <- ifelse(y1 == "yes" & y2 == "no",   1, 0)
  firstLower.n   <- ifelse(y1 == "no"  & y2 == "none", 1, 0)
  # Respondents who face a higher bid in the first question
  firstHigher.y  <- ifelse(y1 == "yes" & y2 == "none", 1, 0)
  firstHigher.ny <- ifelse(y1 == "no"  & y2 == "yes",  1, 0)
  firstHigher.nn <- ifelse(y1 == "no"  & y2 == "no",   1, 0)
} else {
  # Respondents who face a lower bid in the first question
  firstLower.yy  <- ifelse(y1 == 1 & y2 ==  1,  1, 0)
  firstLower.yn  <- ifelse(y1 == 1 & y2 ==  0,   1, 0)
  firstLower.n   <- ifelse(y1 == 0 & y2 == -9, 1, 0)
  # Respondents who face a higher bid in the first question
  firstHigher.y  <- ifelse(y1 == 1 & y2 == -9, 1, 0)
  firstHigher.ny <- ifelse(y1 == 0 & y2 ==  1, 1, 0)
  firstHigher.nn <- ifelse(y1 == 0 & y2 ==  0, 1, 0)
}
# Creating response variables
yy <- firstLower.yy + firstHigher.y
yn <- firstLower.yn + firstHigher.ny
nn <- firstLower.n  + firstHigher.nn
BID <- Llow * (firstLower.yy + firstLower.yn + firstLower.n) +
       Lhig * (firstHigher.y + firstHigher.ny + firstHigher.nn)  # start bid
bid <- cbind(Llow, Lhig)
yvar <- cbind(yy, yn, nn)

#    # Creating a design matrix
# Revised in June 2016
#    bidf <- formula(formula, lhs = 0, rhs = 2)
#    bid <- model.frame(bidf, data)  # the first and the second stage bids
#bid <- model.part(formula, data, lhs = 0, rhs = 2)
#
#    BID <- ifelse(bid[, 1] > bid[, 2], bid[, 2], bid[, 1])
#
#    yvar <- cbind(yy, yn, ny, nn)   # yes/no to "bid"
    
# Revised in June 2016
#    ff2 <- formula(formula, lhs = 0, rhs = 1)
#    X <- model.frame(ff2, data)
#    mmX <- model.matrix(ff2, X)
X <- model.part(formula, data, lhs = 0, rhs = 1)
mmX <- model.matrix(formula, data, lhs = 0, rhs = 1)

    tmp.data <- data.frame(y1, mmX, BID)

   # obtaining initial parameter values by logit model
   if(is.null(par)){
        f.stage <- glm(y1~. -1, family = binomial(link = "probit"), data = tmp.data)
         ini <- f.stage$coefficients # saving initial values for ML estimation
         npar <- length(ini)
         ini[npar] <- ifelse(ini[npar] > 0, -ini[npar], ini[npar])     # gives a negative initial value for the bid coefficient
         if (substr(dist, 1, 4) == "log-" | dist == "weibull") names(ini)[npar] <- "log(bid)"
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
#        ny <- dvar[, 3]
#        nn <- dvar[, 4]
nn <- dvar[, 3]

        X1 <- cbind(ivar, bid[, 1])
        X2 <- cbind(ivar, bid[, 2])
        
        ll <- 
        sum(plogis(-X2[yy == 1, ]%*%param, lower.tail = FALSE, log.p = TRUE))  + 

#        sum(plogis(-X2[nn == 1, ]%*%param, lower.tail = TRUE, log.p = TRUE))   + 

sum(plogis(-X1[nn ==1, ] %*% param, lower.tail = TRUE, log.p = TRUE)) +
        sum(log(plogis(-X2[yn == 1, ]%*%param, lower.tail = TRUE, log.p = FALSE) - 
                plogis(-X1[yn == 1, ]%*%param, lower.tail = TRUE, log.p = FALSE))) #+   

#        sum(log(plogis(-X1[ny == 1, ]%*%param, lower.tail = TRUE, log.p = FALSE) - 
#                 plogis(-X2[ny == 1, ]%*%param, lower.tail = TRUE, log.p = FALSE)))  

        ifelse(is.finite(ll), return(-ll), NaN) 
      }
  } else if(dist == "normal" | dist == "log-normal") {
      # likelihood function
        dbLL <- function(param, dvar, ivar, bid){
          yy <- dvar[, 1]
          yn <- dvar[, 2]
#          ny <- dvar[, 3]
#          nn <- dvar[, 4]
nn <- dvar[, 3]

          X1 <- cbind(ivar, bid[, 1])
          X2 <- cbind(ivar, bid[, 2])
          
          ll <- 
          sum(pnorm(-X2[yy == 1, ]%*%param, lower.tail = FALSE, log.p = TRUE))  + 

#          sum(pnorm(-X2[nn == 1, ]%*%param, lower.tail = TRUE, log.p = TRUE))   + 
sum(pnorm(-X1[nn == 1,] %*% param, lower.tail = TRUE, log.p = TRUE)) +

          sum(log(pnorm(-X2[yn == 1, ]%*%param, lower.tail = TRUE, log.p = FALSE) - 
                   pnorm(-X1[yn == 1, ]%*%param, lower.tail = TRUE, log.p = FALSE))) #+   
#          sum(log(pnorm(-X1[ny == 1, ]%*%param, lower.tail = TRUE, log.p = FALSE) - 
#                   pnorm(-X2[ny == 1, ]%*%param, lower.tail = TRUE, log.p = FALSE)))  
        ifelse(is.finite(ll), return(-ll), NaN) 
        }
  } else if(dist == "weibull"){
        dbLL <- function(param, dvar, ivar, bid){
          yy <- dvar[, 1]
          yn <- dvar[, 2]
#          ny <- dvar[, 3]
#          nn <- dvar[, 4]
nn <- dvar[, 3]
          
          X1 <- cbind(ivar, bid[, 1])
          X2 <- cbind(ivar, bid[, 2])
          
          ll <- 
          sum(pweibull(exp(-X2[yy == 1, , drop = FALSE]%*%param), shape = 1, lower.tail = FALSE, log.p = TRUE))  + 

#          sum(pweibull(exp(-X2[nn == 1, , drop = FALSE]%*%param), shape = 1, lower.tail = TRUE, log.p = TRUE))   + 
sum(pweibull(exp(-X1[nn == 1, , drop = FALSE] %*% param), shape = 1, lower.tail = TRUE, log.p = TRUE))   + 

          sum(log(pweibull(exp(-X2[yn == 1, , drop = FALSE]%*%param), shape = 1, lower.tail = TRUE, log.p = FALSE) - 
                   pweibull(exp(-X1[yn == 1, , drop = FALSE]%*%param), shape = 1, lower.tail = TRUE, log.p = FALSE))) #+   
#          sum(log(pweibull(exp(-X1[ny == 1, , drop = FALSE]%*%param), shape = 1, lower.tail = TRUE, log.p = FALSE) - 
#                   pweibull(exp(-X2[ny == 1, , drop = FALSE]%*%param), shape = 1, lower.tail = TRUE, log.p = FALSE)))  
        ifelse(is.finite(ll), return(-ll), NaN) 
        }
  }
  
  # ML estimation of double-bounded dichotomous choice
  suppressWarnings(
        dbchoice <- optim(ini, fn = dbLL, method="BFGS", hessian = TRUE, dvar = yvar, ivar = mmX, bid = bid, control = list(abstol = 10^(-30)))
    )
  npar <- length(dbchoice$par)
  
  terms <- terms(formula)

# Revised in June 2016
#  fac <- which(attr(attr(X, "terms"), "dataClasses") == "factor")
fac <- which(sapply(X, is.factor) == TRUE)

  xlevels <- as.list(fac)
  j <- 0
  for (i in fac) {
    j <- j + 1
    xlevels[[j]] <- levels(X[[i]])
  }
  contrasts <- attr(mmX, "contrasts")

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
      covariates = mmX,             # a matrix of the covariates
      bid = bid,                    # the suggested bid
      yn = cbind(y1, y2),           # the acceptance/rejection variable
      data.name = data,             # the data matrix
      terms = terms,
      contrasts = contrasts,
# Revised in June 2016
data = original.data,
      xlevels = xlevels)

#  class(output) <- "dbchoice"       # setting the object class
class(output) <- c("oohbchoice", "dbchoice")       # setting the object class

  return(output)
}




##### a function for Kaplan-Meier-Turnbull nonparametric approach to analyze OOHB choice data #####

turnbull.oohb <- function(formula, data, subset, conf.int = FALSE, B = 200, conf.level = 0.95, timeMessage = FALSE, ...){
    if (!inherits(formula, "Formula"))
        formula <- Formula(formula)

    if(missing(data)) data <- environment(formula)
    
    # stop if the LHS does not contain two variables
    if(length(formula[[2]]) != 3) stop("LHS variable in the formula must be like y1 + y2 ")
    
    mf <- match.call(expand.dots = TRUE)
    m <- match(c("formula", "data", "subset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$formula <- formula
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    original.data <- data
    data <- mf

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
    
    nobs <- length(y2)
    
    P1 <- formula[[3]][[2]]   # extracting from the formula the name of the variable for the first bid
    P2 <- formula[[3]][[3]]   # extracting from the formula the name of the variable for the second bid
#    first <- eval(P1, data)       # the first bids
#    second <- eval(P2, data)      # the second bids
Llow <- eval(P1, data)
Lhig <- eval(P2, data)
    
    # making dummy variables for the first and second bids
#    if(is.factor(y1)){   # when the yes/no variables are defined as factor
#        yy <- ifelse(y1 == "yes" & y2 == "yes", 1, 0)
#        yn <- ifelse(y1 == "yes" & y2 == "no", 1, 0)
#        ny <- ifelse(y1 == "no" & y2 == "yes", 1, 0)
#        nn <- ifelse(y1 == "no" & y2 == "no", 1, 0)
#    } else {
#        yy <- ifelse(y1 == 1 & y2 == 1, 1, 0)
#        yn <- ifelse(y1 == 1 & y2 == 0, 1, 0)
#        ny <- ifelse(y1 == 0 & y2 == 1, 1, 0)
#        nn <- ifelse(y1 == 0 & y2 == 0, 1, 0)
#    }
if (is.factor(y1)) {
  # Respondents who face a lower bid in the first question
  firstLower.yy  <- ifelse(y1 == "yes" & y2 == "yes",  1, 0)
  firstLower.yn  <- ifelse(y1 == "yes" & y2 == "no",   1, 0)
  firstLower.n   <- ifelse(y1 == "no"  & y2 == "none", 1, 0)
  # Respondents who face a higher bid in the first question
  firstHigher.y  <- ifelse(y1 == "yes" & y2 == "none", 1, 0)
  firstHigher.ny <- ifelse(y1 == "no"  & y2 == "yes",  1, 0)
  firstHigher.nn <- ifelse(y1 == "no"  & y2 == "no",   1, 0)
} else {
  # Respondents who face a lower bid in the first question
  firstLower.yy  <- ifelse(y1 == 1 & y2 ==  1,  1, 0)
  firstLower.yn  <- ifelse(y1 == 1 & y2 ==  0,   1, 0)
  firstLower.n   <- ifelse(y1 == 0 & y2 == -9, 1, 0)
  # Respondents who face a higher bid in the first question
  firstHigher.y  <- ifelse(y1 == 1 & y2 == -9, 1, 0)
  firstHigher.ny <- ifelse(y1 == 0 & y2 ==  1, 1, 0)
  firstHigher.nn <- ifelse(y1 == 0 & y2 ==  0, 1, 0)
}

#    left <- ifelse(yy == 1 | ny == 1, second, ifelse(yn == 1, first, 0))    # lower bound of WTP
#    right <- ifelse(yn ==1 | nn == 1, second, ifelse(ny == 1, first, Inf))  # upper bound of WTP

left  <- ifelse(firstLower.yy == 1 | firstHigher.y  == 1, Lhig,
         ifelse(firstLower.yn == 1 | firstHigher.ny == 1, Llow, 0))
right <- ifelse(firstLower.yn == 1 | firstHigher.ny == 1, Lhig,
         ifelse(firstLower.n  == 1 | firstHigher.nn == 1, Llow, Inf))

    unq.bid <- sort(unique(c(left, right)))   # unique bids including Inf
    
    # estimating nonparametric survival function. icfit function is defined in interval package
    turnbull <- icfit(L = left, R = right, conf.int = conf.int, 
                      control = icfitControl(timeMessage = timeMessage, B = B,
                                             conf.level = conf.level))
    
    # arranging outcomes into a single list variable
    output <- list(
        left = left,
        right = right,
        turnbull = turnbull,
        unq.bid = unq.bid
    )
    
    class(output) <- "turnbull"
    return(output)
    
}

######################################################################
