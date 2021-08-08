sbspike <- function(
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

  if (!inherits(formula, "formula")) {
    stop("invalid formula")
  }

  if(missing(data)) {
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
  mm.data <- model.matrix(formula, data = data, rhs = 1:2)

# defining a dependent variable 
  y1    <- model.part(formula, data = data, lhs = 1)[[1]]
  spk   <- model.part(formula, data = data, lhs = 2)[[1]]

  y1spk <- 1 * (y1 == 1) + 2 * (y1 == 0) *(spk == 1) + 3 * (spk == 0)

# creating a design matrix
  BID <-   model.part(formula, data, lhs = 0, rhs = 2)

  X   <-   model.part(formula, data, lhs = 0, rhs = 1)
  mmX <- model.matrix(formula, data, lhs = 0, rhs = 1)

  if(!any(colnames(mmX) == "(Intercept)")) {
    stop(message = "constant (intercept) term is required for the formula")
  }

# obtaining initial parameter values
  tmp.data <- data.frame(y1, mm.data[, -1, drop = FALSE])

  if(is.null(par)) {
    f.stage <- glm(y1 ~. , family = binomial(link = "logit"),
                   data = tmp.data)
    ini.par <- f.stage$coefficients
  } else {
    if(length(par) != ncol(tmp.data)) {
      stop("length of 'par' must be equal to number of independent variables")
    }
    ini.par <- par
    f.stage <- ini.par
  }

# defining a log-likelihood function
  sbLL <- function(param, dvar, ivar) {
    y1spk <- dvar  # 1 = yes; 2 = no-yes; 3 = no-no
    X0    <- X1 <- ivar
    X0[, ncol(X0)] <- 0 # Bid variable in X0 has the value of 0

    ll <- sum(log(plogis(-X1[y1spk == 1, , drop = FALSE] %*% param,
                         lower.tail = FALSE, log.p = FALSE))) + 
          sum(log(plogis(-X1[y1spk == 2, , drop = FALSE] %*% param,
                         lower.tail = TRUE,  log.p = FALSE) -
                  plogis(-X0[y1spk == 2, , drop = FALSE] %*% param,
                         lower.tail = TRUE,  log.p = FALSE))) +
          sum(log(plogis(-X0[y1spk == 3, , drop = FALSE] %*% param,
                         lower.tail = TRUE,  log.p = FALSE)))

    ifelse(is.finite(ll), return(-ll), NaN)
  }

# estimating SBDC spike model
  suppressWarnings(
    optim.out <- optim(ini.par,
                       fn = sbLL,
                       method = "BFGS",
                       hessian = TRUE,
                       dvar = y1spk,
                       ivar = mm.data))

# storing factor levels
  fac <- which(sapply(X, is.factor) == TRUE)
  xlevels <- as.list(fac)
  j <- 0
  for (i in fac) {
    j <- j + 1
    xlevels[[j]] <- levels(X[[i]])
  }

# storing outcomes into a list object
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
    bid          = BID,
    yn           = y1,
    data.name    = data,
    terms        = terms(formula),
    contrasts    = attr(mmX, "contrasts"),
    data         = original.data,
    xlevels      = xlevels)

# setting the object class
  class(output) <- c("spike", "sbspike")

  return(output)
}

