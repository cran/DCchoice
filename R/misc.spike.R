###############################################################################
print.spike <- function(
  x,
  digits = max(3, getOption("digits") - 1),
  ...)
{
  if(!x$convergence) {
    cat("The optimization did not converge\n")
  }
  cat("\nSpike model\n")
  print.default(format(x$coef, digits = digits), print.gap = 1, quote = FALSE)

  invisible(x)
}

###############################################################################
summary.spike <- function(
  object,
  ...)
{
  # obtaining necessary components from the object
  coef <- object$coefficients
  npar <- length(coef)
  se <- sqrt(diag(solve(object$Hessian)))
  bid <- object$bid
  X <- object$covariates

  # estimating the null model
  formula_null <- formula(object$formula)
  formula_null[[3]][[2]] <- 1
  if(inherits(object, "sbspike")){
    null_model <-   sbspike(formula_null,
                            data = eval(object$data),
                            par = coef[c(1, npar)])
  } else if(inherits(object, "dbspike")) {
    null_model <-   dbspike(formula_null,
                            data = eval(object$data),
                            par = coef[c(1, npar)])
  } else if(inherits(object, "oohbspike")) {
    null_model <- oohbspike(formula_null,
                            data = eval(object$data),
                            par = coef[c(1, npar)])
  } else {
    stop("object should be S3 class 'spike'")
  }

  # function for obrtaining AIC and BIC
  akaike <- function(loglik, npar, k ){
    -2*loglik + k*npar
  }

  # computing WTPs
  WTPs <- wtp_spike(object = object)
  object$medianWTP <- WTPs$medianWTP
  object$meanWTP <- WTPs$meanWTP
  
  # conducting likelihood ratio test
  LR <- -2*(null_model$loglik - object$loglik)
  d.f <- length(object$coefficients) - length(null_model$coefficients)
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

  class(object) <- "summary.spike"

  return(object)
}

###############################################################################
print.summary.spike <- function(
  x,
  digits = max(3, getOption("digits") - 1),
  ...)
{
  cat("\nCall:",
      deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)),
      "", sep = "\n")
  cat("Formula:",
      deparse(x$formula, width.cutoff = floor(getOption("width") * 0.85)),
      "", sep = "\n")

  if(!x$convergence) {
    cat("The optimization did not converge\n")
  }
  cat("Coefficients:", "\n")
  printCoefmat(x$coef, digits = 4, dig.tst = 4)
  
  cat("\nDistribution: spike logistic")
  cat("\nNumber of Obs.:", formatC(x$nobs, digits = 0), "\n")
  cat("Log-likelihood:", formatC(x$loglik, format="f", digits = digits), "\n")
  cat("\nLR statistic:", formatC(x$LR.test[1], format="f", digits = 3), 
      "on", formatC(x$LR.test[2], digits = 0), 
      "DF, p-value:", formatC(x$LR.test[3], format="f", digits = 3), "\n")
  cat("AIC:", formatC(x$AIC[1], format="f", digits = digits), 
      ", BIC:", formatC(x$AIC[2], format="f", digits = digits), "\n")
  cat("\nIterations:", formatC(x$niter, digits = 0), "")
  cat("\nConvergence:", x$convergence, "\n")

  cat("\nWTP estimates:")
  if(is.finite(x$meanWTP)){
    cat("\n Mean :", x$meanWTP, "", sep = " ")
  } else {
    cat("\n Mean :", x$meanWTP, "(because of |beta_Lbid| < 1)", sep = " ")
  }
  cat("\n Median:", x$medianWTP, "\n", sep = " ")

}

###############################################################################
spikeCoef <- function(
  object)
{
  X <- object$covariates
  B <- object$coefficients

  names(B) <- NULL
  ncoefs   <- length(B) 
  mX       <- colMeans(X)
  XB       <- sum(mX * B[-ncoefs])
  VCOV     <- vcov(object)[-ncoefs, -ncoefs]
  dlogisXB <- dlogis(XB)

  spike  <- 1/(1 + exp(XB))
  spikese <- as.vector(sqrt(dlogisXB %*% t(mX) %*% VCOV %*% mX %*% dlogisXB))

  # creating a table for coefficients, se, etc. 
  z <- spike/spikese
  p <- round(2 * pnorm(-abs(z)), 6)
  rtn <- cbind(spike, spikese, z, p)
  colnames(rtn) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  rownames(rtn) <- "Spike"

  rtn
}

###############################################################################
wtp_spike <- function(
  object,
  b = NULL) 
{
  if(inherits(object, "spike")) {
    X <- object$covariates
    B <- object$coefficients
  } else {
    X <- object
    B <- b
  }

  names(B) <- NULL
  ncoefs   <- length(B) 
  Bbid     <- B[ncoefs]
  XB <- sum(colMeans(X) * B[-ncoefs])

  if ((1/(1 + exp(-XB))) < 0.5) {
    medianWTP <- 0
  } else {
    medianWTP <- -XB/Bbid
  }

  if (Bbid > 0) {
    meanWTP <- NA
  } else {
    meanWTP <- -log(1+exp(XB))/Bbid
  }

  list(meanWTP = meanWTP, medianWTP = medianWTP)
}

###############################################################################
logLik.spike <- function(
  object,
  ...) 
{
  rtn <- object$loglik
  rtn
}

###############################################################################
vcov.spike <- function(
  object,
  ...) 
{
  solve(object$Hessian)
}

###############################################################################
plot.spike <- function(
  x,
  type = "l",
  main = NULL,
  sub = NULL, 
  xlab = "Bid",
  ylab = "Probability",
  lwd = 3,
  lty = 1,
  xlim = c(0, max(x$bid)), 
  ylim = c(0, 1),
  bid = c(0, max(x$bid)),
  ...) 
{
  minbid <- min(bid)
  maxbid <- max(bid)
  bidseq <- seq(minbid, maxbid, by = (maxbid - minbid)/100)
  b      <- x$coefficients
  npar   <- length(b)
  Xcov   <- colMeans(x$covariates)
  Vcov   <- sum(Xcov * b[-npar])
  V      <- Vcov + bidseq * b[npar]
  pr     <- plogis(-V, lower.tail = FALSE, log.p = FALSE)

  plot(x = c(0, bidseq),
       y = c(1, pr),
       xlab = xlab,
       ylab = ylab,
       main = main,
       lwd = lwd,
       type = type,
       xlim = xlim,
       ylim = ylim,
       axes = FALSE,
       ...) 
  axis(side = 1, pos = 0)
  axis(side = 2, pos = 0)

  invisible(list(x = bidseq, utility = V, probability = pr))
}

###############################################################################
ct2df <- function(
  x,
  bid1  = "bid1",  
  bid2h = "bidh",  
  bid2l = "bidl",
  bidl  = "bidl",
  bidh  = "bidh",
  nny   = "nny",
  nnn   = "nnn",  
  yy    = "yy",
  yn    = "yn",
  ny    = "ny",
  nn    = "nn",
  y     = "y", 
  n     = "n",
  n_y   = "n_y",
  n_n   = "n_n",
  nn_y  = "nn_y",
  nn_n  = "nn_n",
  type  = "double",
  spike = FALSE) 
{
  if(isTRUE(spike)) {
    ct2df.spike(
      x     = x,
      bid1  = bid1,  
      bid2h = bid2h,  
      bid2l = bid2l,
      bidl  = bidl,
      bidh  = bidh,
      nny   = nny,
      nnn   = nnn,  
      yy    = yy,
      yn    = yn,
      ny    = ny,
      nn    = nn,
      y     = y, 
      n     = n,
      n_y   = n_y,
      n_n   = n_n,
      nn_y  = nn_y,
      nn_n  = nn_n,
      type  = type)
  } else {
    ct2df.default(
      x     = x,
      bid1  = bid1,
      bid2h = bid2h,
      bid2l = bid2l, 
      yy    = yy,
      yn    = yn, 
      ny    = ny,
      nn    = nn, 
      y     = y,
      n     = n,
      type  = type)
  }
}

###############################################################################
ct2df.spike <- function(
  x,
  bid1,  
  bid2h,  
  bid2l,
  bidl,
  bidh,
  nny,
  nnn,  
  yy,
  yn,
  ny,
  nn,
  y, 
  n,
  n_y,
  n_n,
  nn_y,
  nn_n,
  type) {
  if (type == "single") {
    reshape.x <- x[c(bid1, y, ny, nn)]

    bid1 <- c(rep(reshape.x[, 1], reshape.x[, 2]),
              rep(reshape.x[, 1], reshape.x[, 3]),
              rep(reshape.x[, 1], reshape.x[, 4]))

    R1   <- c(rep(1L, sum(reshape.x[, 2])),
              rep(0L, sum(reshape.x[, 3])),
              rep(0L, sum(reshape.x[, 4])))

    S    <- c(rep(1L, sum(reshape.x[, 2])),
              rep(1L, sum(reshape.x[, 3])),
              rep(0L, sum(reshape.x[, 4])))

    cv.data <- data.frame(R1 = R1, S = S, bid1 = bid1)
  } else if (type == "double") {
  # double-bounded
    reshape.x <- x[c(bid1, bid2h, bid2l, yy, yn, ny, nny, nnn)] 
    colnames(reshape.x) <- 
      c("B1", "B2H", "B2L", "yy", "yn", "ny", "nny", "nnn")

    bid.table      <- reshape.x[, c(1, 2, 3)]
    response.table <- reshape.x[, c(4, 5, 6, 7, 8)]

    B1 <- rep(bid.table[, 1], rowSums(response.table))
    R  <- rep(names(response.table), response.table[1, ])
    for (i in 2:nrow(x)) {
      R <- c(R, rep(names(response.table), response.table[i, ]))
    }

    data <- data.frame(
      B1 = B1,
      R  = R,
      R1 = c(R == "yy") + c(R == "yn"),
      R2 = c(R == "yy") + c(R == "ny"),
      S  = c(R != "nnn")*1)

    cv.data <- merge(bid.table, data,  by = "B1")
    cv.data$bid1 <- cv.data$B1
    cv.data$bid2 <- cv.data$B2H * cv.data$R1 + 
                    cv.data$B2L * (cv.data$R1 == 0)

  } else {
    reshape.x <- x[c(bidl, bidh, yy, yn, n_y, n_n, y, ny, nn_y, nn_n)]
    colnames(reshape.x) <- 
      c("BL", "BH", "yy", "yn", "n_y", "n_n", "y", "ny", "nn_y", "nn_n")

    bid.table      <- reshape.x[, c(1, 2)]
    response.table <- reshape.x[, c(3, 4, 5, 6, 7, 8, 9, 10)]

    BL <- rep(bid.table[, 1], rowSums(response.table))
    R  <- rep(names(response.table), response.table[1, ])
    for (i in 2:nrow(x)) {
      R <- c(R, rep(names(response.table), response.table[i, ]))
    }

    data <- data.frame(
      BL = BL,
      R  = R,
      R1 = (R == "yy") + (R == "yn") + (R == "y"),
      R2 = (R == "yy") + (R == "ny") - 
           9 * ((R == "n_y") + (R == "n_n") + (R == "y")),
      S  = (R == "yy") + (R == "yn") + (R == "n_y") + 
           (R == "y") + (R == "ny") + (R == "nn_y"))

    cv.data <- merge(bid.table, data,  by = "BL")
  }

  return(cv.data)
}

###############################################################################
bootCI <- function(
  obj,
  nboot = 1000,
  CI = 0.95,
  individual = NULL) 
{
  if(inherits(obj, "spike")) {
      bootCI.spike(obj = obj, nboot = nboot, CI = CI, individual = individual)
  } else {
    bootCI.default(obj = obj, nboot = nboot, CI = CI, individual = individual)
  }
}

###############################################################################
krCI <- function(
  obj,
  nsim = 1000,
  CI = 0.95,
  individual = NULL) 
{
  if(inherits(obj, "spike")) {
      krCI.spike(obj = obj, nsim = nsim, CI = CI, individual = individual)
  } else {
    krCI.default(obj = obj, nsim = nsim, CI = CI, individual = individual)
  }
}

###############################################################################
bootCI.spike <- function (
  obj,
  nboot,
  CI,
  individual)
{
# checking arguments
  if(!inherits(obj, c("spike"))){
    stop("the object must be spike class")
  }
  if(CI > 1 | CI < 0) {
    stop("CI must be between 0 and 1")
  }

# setting objects
  tmp.dat     <- eval(obj$data, parent.frame())
  nobs        <- obj$nobs
  ind         <- 1:nobs
  boot.mean   <- numeric(nboot)
  boot.median <- numeric(nboot)
  fr          <- obj$formula

  if (!is.null(individual)) {
    formula <- formula(obj$formula, lhs = 0, rhs = -2)
    mf.nexX <- model.frame(formula, individual, xlev = obj$xlevels)
    mm.newX <- model.matrix(formula, mf.nexX, contrasts.arg = obj$contrasts)
  }

# simulating and storing mean and median WTPs
  if(inherits(obj, c("dbspike", "oohbspike"))){
    for(i in 1:nboot){
      ind.boot <- sample(ind, nobs, replace = TRUE)
      boot.dat <- tmp.dat[ind.boot, ]  # resampling data
  
      if(!inherits(obj, "oohbspike")){
        suppressWarnings(
          tmp.re <-   dbspike(fr, data = boot.dat, par=obj$coefficients))
      } else if(inherits(obj, "oohbspike")) {
        suppressWarnings(
          tmp.re <- oohbspike(fr, data = boot.dat, par=obj$coefficients))
      }
  
      if(tmp.re$convergence){
        if (is.null(individual)) {
          boot <- wtp_spike(object = tmp.re$covariates,
                           b      = tmp.re$coefficients)
        } else {
          boot <- wtp_spike(object = mm.newX,
                           b      = tmp.re$coefficients)
        }
        boot.mean[i]   <- boot$meanWTP
        boot.median[i] <- boot$medianWTP
      } else {
        i < i - 1
      }
    }
  } else if(inherits(obj, "sbspike")){
    for(i in 1:nboot){
      ind.boot <- sample(ind, nobs, replace = TRUE)
      boot.dat <- tmp.dat[ind.boot, ]
      suppressWarnings(
        tmp.re <- sbspike(fr, data = boot.dat, par=obj$coefficients))
      if(tmp.re$convergence){
        if (is.null(individual)) {
          boot <- wtp_spike(object = tmp.re$covariates,
                           b      = tmp.re$coefficients)
        } else {
          boot <- wtp_spike(object = mm.newX,
                           b      = tmp.re$coefficients)
        }
        boot.mean[i]   <- boot$meanWTP
        boot.median[i] <- boot$medianWTP
      } else {
        i < i -1
      }
    }
  }
  
  output <- list(mWTP = boot.mean, medWTP = boot.median)
  
# sorting the simulation outcomes
  boot.mean   <- sort(boot.mean)
  boot.median <- sort(boot.median)

  lb <- 0.5*(1 - CI)
  ub <- CI + lb

  int <- c(ceiling(nboot*lb), floor(nboot*ub))

# 100*CI% bootstrap CI
  CImat <- rbind(boot.mean[int], boot.median[int]) 
  
# the mean estimates in the original outcome
  if (is.null(individual)) {
    tmp.sum <- wtp_spike(object = obj$covariates,
                        b      = obj$coefficients)
  } else {
    tmp.sum <- wtp_spike(object = mm.newX,
                        b      = obj$coefficients)
  }
  
  out <- cbind(c(tmp.sum$meanWTP, tmp.sum$medianWTP), CImat)
  rownames(out) <- c("Mean", "Median")
  colnames(out) <- c("Estimate", "LB", "UB")

  output$out <- out
  
  class(output) <- "bootCI"
  return(output)
}

###############################################################################
krCI.spike <- function(
  obj,
  nsim,
  CI,
  individual)
{
# checking arguments
  if(CI > 1) {
    stop("CI must be between 0 and 1")
  }
  if(!inherits(obj, "spike")) {
    stop("the object must be spike class")
  }

# setting objects  
  X       <- obj$covariates             # the covariates
  npar    <- length(obj$coefficients)   # the number of coefficients
  coef    <- obj$coefficients           # coefficient estimates
  formula <- formula(obj$formula, lhs = 0, rhs = -2)
  VCOV    <- solve(obj$Hessian)

# simulating coefficients
  kr.coef    <- t(mvrnorm(nsim, coef, VCOV))   # mvrnorm in MASS package
  kr.coefbid <- kr.coef[npar,  , drop = FALSE] # coefficient on the bid
  kr.coefcov <- kr.coef[-npar, , drop = FALSE] # coefficients on the covariates

  if(is.null(individual)) {
    kr.Xb <- colSums(colMeans(X[, -npar, drop = FALSE]) * kr.coefcov)
  } else {
    mf.nexX <- model.frame(formula, individual, xlev = obj$xlevels)
    mm.newX <- model.matrix(formula, mf.nexX, contrasts.arg = obj$contrasts)
    newX    <- as.vector(mm.newX)
    kr.Xb   <- colSums(newX * kr.coefcov)
  }

# calculating and storing mean and median WTPs
  kr.median  <- ifelse((1 + exp(-kr.Xb))^-1 < 0.5, 0, -kr.Xb/kr.coefbid)
  kr.meanWTP <- ifelse(kr.coefbid < 0, -log(1 + exp(kr.Xb))/kr.coefbid, NA)
  output <- list(mWTP = kr.meanWTP, medWTP = kr.median)
  
# sorting the simulation outcomes
  kr.meanWTP <- sort(kr.meanWTP, na.last = TRUE)
  kr.median  <- sort(kr.median)

  lb <- 0.5*(1 - CI)   # lower bound of the empirical distribution
  ub <- CI + lb        # upper bound of the empirical distribution
  int <- c(ceiling(nsim*lb), floor(nsim*ub))

# 100*CI% simulated conficdence intervals
  CImat <- rbind(kr.meanWTP[int], kr.median[int])

# the mean estimates in the original outcome
  if (is.null(individual)) {
    tmp.sum <- wtp_spike(object = obj$covariates, b = obj$coefficients)
  } else {
    tmp.sum <- wtp_spike(object = mm.newX, b = obj$coefficients)
  }
 
  out <- cbind(c(tmp.sum$meanWTP, tmp.sum$medianWTP), CImat)
  rownames(out) <- c("Mean", "Median")
  colnames(out) <- c("Estimate", "LB", "UB")

  output$out <- out
  
  class(output) <- "krCI"
  return(output)
}

###############################################################################

