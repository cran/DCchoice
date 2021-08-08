oohbspike <- function (
  formula,
  data,
  subset,
  na.action = na.omit, 
  par = NULL,
  ...)
{
# storing a function call
  cl <- match.call()

# checking arguments
  if (!inherits(formula, "Formula")){
    formula <- Formula(formula)
  }

  if (!inherits(formula, "formula")){
    stop("invalid formula")
  }

  if(length(formula[[2]]) != 3){
    stop("LHS variable in the formula must be like y1 + y2 ")
  }

  if(missing(data)){
    stop("data frame must be assigned to argument 'data'")
  }

# preparing a data set for analysis
  mf <- match.call(expand.dots = TRUE)
  m  <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$formula <- formula
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  original.data <- data
  data <- mf
  
# defining a dependent variable
  y1  <- model.part(formula, data = data, lhs = 1)[[1]]
  y2  <- model.part(formula, data = data, lhs = 1)[[2]]
  spk <- model.part(formula, data = data, lhs = 2)[[1]]

  Llow <- model.matrix(formula, data = data, rhs = 2)[, 2]
  Lhig <- model.matrix(formula, data = data, rhs = 2)[, 3]

  if (is.factor(y1)) {
  # respondents who face a lower bid in the first question
    firstLower.yy  <- ifelse(y1 == "yes" & y2 == "yes",  1, 0)
    firstLower.yn  <- ifelse(y1 == "yes" & y2 == "no",   1, 0)
    firstLower.n   <- ifelse(y1 == "no"  & y2 == "none", 1, 0)
  # respondents who face a higher bid in the first question
    firstHigher.y  <- ifelse(y1 == "yes" & y2 == "none", 1, 0)
    firstHigher.ny <- ifelse(y1 == "no"  & y2 == "yes",  1, 0)
    firstHigher.nn <- ifelse(y1 == "no"  & y2 == "no",   1, 0)
  } else {
  # respondents who face a lower bid in the first question
    firstLower.yy  <- ifelse(y1 == 1 & y2 ==  1, 1, 0)
    firstLower.yn  <- ifelse(y1 == 1 & y2 ==  0, 1, 0)
    firstLower.n   <- ifelse(y1 == 0 & y2 == -9, 1, 0)
  # respondents who face a higher bid in the first question
    firstHigher.y  <- ifelse(y1 == 1 & y2 == -9, 1, 0)
    firstHigher.ny <- ifelse(y1 == 0 & y2 ==  1, 1, 0)
    firstHigher.nn <- ifelse(y1 == 0 & y2 ==  0, 1, 0)
  }

  yy <- firstLower.yy + firstHigher.y
  yn <- firstLower.yn + firstHigher.ny
  nn <- firstLower.n  + firstHigher.nn

  yvar <- cbind(yy, yn, nn, spk)

# creating a design matrix
  BID <- Llow * (firstLower.yy + firstLower.yn + firstLower.n) +
         Lhig * (firstHigher.y + firstHigher.ny + firstHigher.nn)
  bid <- cbind(Llow, Lhig)

  X   <-   model.part(formula, data, lhs = 0, rhs = 1)
  mmX <- model.matrix(formula, data, lhs = 0, rhs = 1)

  if(!any(colnames(mmX) == "(Intercept)")) {
    stop(message = "constant (intercept) term is required for the formula")
  }

# obtaining initial parameter values
  tmp.data <- data.frame(y1, mmX, BID)

  if(is.null(par)){
    f.stage <- glm(y1~. -1, family = binomial(link = "logit"), 
                   data = tmp.data)
    ini     <- f.stage$coefficients
    npar    <- length(ini)
    ini[npar] <- ifelse(ini[npar] > 0, -ini[npar], ini[npar])
    names(ini)[1] <- "(Intercept)"
  } else {
    if(length(par) != ncol(tmp.data)-1) {
      stop("length of 'par' must be equal to number of independent variables")
    }
    ini <- par
    f.stage <- ini
  }

# defining a log-likelihood function
  oohbLL <- function(param, dvar, ivar, bid){
    yy  <- dvar[, 1]
    yn  <- dvar[, 2]
    nn  <- dvar[, 3]
    spk <- dvar[, 4]
    nny <- nn * (spk == 1)
    nnn <- nn * (spk == 0)
    X1  <- cbind(ivar, bid[, 1])
    X2  <- cbind(ivar, bid[, 2])
    X0  <- X1
    X0[, ncol(X0)] <- 0

    ll <- sum(    plogis(-X2[yy  == 1, , drop = FALSE] %*% param,
                         lower.tail = FALSE, log.p =  TRUE))  + 
          sum(log(plogis(-X2[yn  == 1, , drop = FALSE] %*% param,
                         lower.tail = TRUE,  log.p = FALSE) - 
                  plogis(-X1[yn  == 1, , drop = FALSE] %*% param,
                         lower.tail = TRUE,  log.p = FALSE))) +   
          sum(log(plogis(-X1[nny == 1, , drop = FALSE] %*% param,
                         lower.tail = TRUE,  log.p = FALSE) -
                  plogis(-X0[nny == 1, , drop = FALSE] %*% param,
                         lower.tail = TRUE,  log.p = FALSE))) +
          sum(    plogis(-X0[nnn == 1, , drop = FALSE] %*% param,
                         lower.tail = TRUE,  log.p =  TRUE))

    ifelse(is.finite(ll), return(-ll), NaN) 
  }
  
# estimating OOHB spike model
  suppressWarnings(
    optim.out <- optim(ini,
                       fn = oohbLL,
                       method = "BFGS",
                       hessian = TRUE,
                       dvar = yvar,
                       ivar = mmX,
                       bid = bid,
                       control = list(abstol = 10^(-30))))

# storing factor levels
  fac <- which(sapply(X, is.factor) == TRUE)
  xlevels <- as.list(fac)
  j <- 0
  for (i in fac) {
    j <- j + 1
    xlevels[[j]] <- levels(X[[i]])
  }

# arranging outcomes into a list object
  output <- list(
    f.stage      = f.stage,
    optim.out    = optim.out,
    coefficients = optim.out$par,
    call         = cl,
    formula      = formula,
    Hessian      = optim.out$hessian,
    loglik       = -optim.out$value,
    convergence  = ifelse(optim.out$convergence == 0, TRUE, FALSE),
    niter        = optim.out$counts,
    nobs         = length(y1),
    covariates   = mmX,
    bid          = bid,
    yn           = cbind(y1, y2),
    data.name    = data,
    terms        = terms(formula),
    contrasts    = attr(mmX, "contrasts"),
    data         = original.data,
    xlevels      = xlevels)

# setting the object class
  class(output) <- c("spike", "oohbspike")

  return(output)
}

