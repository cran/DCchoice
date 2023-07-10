dccvResponse <- function(bid, b, n, seed = NULL) {
  if(is.vector(bid)) {
    nBids <- length(bid)
  } else if(is.matrix(bid)) {
    nBids <- nrow(bid)
    if(ncol(bid) != 3)
      stop("bid in matrix should consist of three columns")
  } else {
    stop("bid should be vector or matrix")
  }

  if(n %% nBids != 0)
   if(is.vector(bid)) {
     stop("n should be divisible by number of bids")
   } else {
     stop("n should be divisible by number of rows of bid matrix")
   }

  nPerBid <- n/nBids

  if(is.vector(bid)){
    designMat <- cbind(1, rep(bid, each = nPerBid))
    set.seed(seed)
    e <- rlogis(n)
    Y <- as.integer(designMat %*% b + e > 0)
    rtn <- data.frame(bid1 = designMat[, -1], R1 = Y)
  } else {
    designMat <- cbind(1, bid[rep(x = seq_len(nBids), each = nPerBid), ])
    set.seed(seed)
    e <- rlogis(n)
    Y1  <- as.integer(designMat[, c(1, 2)] %*% b + e > 0)
    Y2H <- as.integer(designMat[, c(1, 3)] %*% b + e > 0)
    Y2L <- as.integer(designMat[, c(1, 4)] %*% b + e > 0)
    Y2H[Y1 == 0] <- 0
    Y2L[Y1 == 1] <- 0
    rtn <- data.frame(
             bid1 = designMat[, 2],
             bid2 = designMat[, 3] * (Y1 == 1) + designMat[, 4] * (Y1 == 0),
             R1   = Y1,
             R2   = Y2H + Y2L)
  }

  rtn
}

