######################################################################################
## mPLPnp: an R package for efficient inference in nonparametric regression
## Authors: Giuseppe Cavaliere, Sílvia Gonçalves, Morten Ø. Nielsen, & Edoardo Zanelli
## version: 1.0 (10.12.2025)
######################################################################################


mPLPnp <- function(y, x, h, eval = NULL, p = NULL, kersel = NULL, res = NULL, alpha = NULL, fast = NULL, bnd = NULL, g.loo = NULL) {
  
  # y       <- dependent variable
  # x       <- regressor
  # h       <- bandwidth
  # eval    <- evaluation point 
  # p       <- polynomial order                 (default = 1)
  # kersel  <- kernel function used             (default = "tri")
  # res     <- residuals class used             (default = "cct-hc3")
  # alpha   <- significance level               (default = 0.05)
  # fast    <- fast std. error evaluation       (default = FALSE)
  # bnd     <- nature of the eval (int or bnd)  (default = FALSE)
  # g.loo   <- leave-one-out Bias estimation    (default = FALSE)
  
  
  ############################# Error Checks ############################# 
  
  if (!is.numeric(h) || h <= 0) {
    stop("'h' must be a positive number")
  }
  
  
  if (!is.numeric(eval) & !is.null(eval)) {
    stop("'eval' must be numeric")
  } else if (is.null(eval)) {
    eval <- 0
  }
  
  
  if (length(p) == 0) {
    p  <- 1
  } else if (length(p) > 1) {
    stop("Polynomial order p incorrectly specified.\n")
  } 
  
  
  if (!is.null(kersel) && is.character(kersel) && length(kersel) == 1 && !(kersel %in% c("uni", "uniform", "tri", "triangular", "epa", "epanechnikov"))){
    stop("Kernel function incorrectly specified.\n")
  } else if (is.null(kersel)){
    kersel <- "tri"
  }
  
  if (!is.null(res) && is.character(res) && length(res) == 1 && !(res %in% c("loo","cct-hc0","cct-hc1","cct-hc2","cct-hc3"))) {
    stop("Residuals incorrectly specified.\n")
  } else if (is.null(res)) {
    res <- "cct-hc3"
  }
  
  
  if (is.null(alpha))  {
    alpha  <- 0.05
  } else if (alpha <= 0 | alpha >=1) {
    stop("significance level incorrectly specified.\n")
  }
  
  
  if (is.null(fast)) {
    fast <- FALSE
  } else if (!is.logical(fast) || length(fast) > 1) {
    stop("'fast' must be a single TRUE/FALSE")
  }
  
  if (is.null(bnd)) {
    bnd <- FALSE
  } else if (!is.logical(bnd) || length(g.loo) > 1) {
    stop("'bnd' must be a single TRUE/FALSE")
  }
  
  
  if (is.null(g.loo)) {
    g.loo <- FALSE
  } else if (!is.logical(g.loo) || length(g.loo) > 1) {
    stop("'g.loo' must be a single TRUE/FALSE")
  }
  
  
  ##########################  Initialization  ########################## 
  
  # Sample size
  n <- length(x)
  
  
  # Bw regularization   
  bwcheck <- 21
  if (!is.null(bwcheck)) {
    bw.min   <- sort(abs(x-eval))[bwcheck]
    h        <- max(h, bw.min)
  }
  
  
  ##########################  function and Bias estimation    ########################## 
  
  # ATE estimation
  ghat   <- LP.est(y=y, x=x, x0=eval, h=h, p=p, v=0, kersel=kersel)
  
  # Initialization of bias estimation
  if (res == "loo") {U <- (x-eval)/(2*h)} else {U <- (x-eval)/(h)}    
  U.Z              <- abs(U) < 1 
  nonzero_indices  <- which(U.Z != 0) 
  ghat_vec <- numeric(n)
  if (g.loo == TRUE) ghat_vec.loo <- numeric(n)
  if (fast == FALSE) C.xi <- numeric(n)
  
  # Loop for bias estimation
  for (ix in seq_along(nonzero_indices)) {
    i <- nonzero_indices[ix]
    x0i <- x[i]
    ghat_vec[i]     <- LP.est(y=y,x=x,x0=x0i,h=h,p=p,v=0,kersel=kersel)
      
    if (fast == FALSE) C.xi[i]  <- LP.est(y=((x-x0i)/h)^(p+1), x=x, x0=x0i, h=h, p=p, v=0, kersel=kersel)
    if (g.loo == TRUE) {
      x.loo <- x[-i];   yp.loo <- y[-i]
      ghat_vec.loo[i] <- LP.est(y=y.loo, x=x.loo, x0=x0i, h=h, p=p, v=0, kersel=kersel)
    }
  }
  
  # Estimation of Q at each side of the cutoff
  if (fast == FALSE) {
    C.LP <- LP.est(y=C.xi,x=x,x0=eval,h=h,p=p,v=0,kersel=kersel)
    C    <- LP.est(y=((x-eval)/h)^(p+1),x=x,x0=eval,h=h,p=p,v=0,kersel=kersel)
    Q    <- C/C.LP;  
  } else if (fast == TRUE) {
    if (bnd == TRUE) {
      if (kersel == "tri") {
        Q <- 1.4082
      } else if (kersel == "epa") {
        Q <- 1.3571
      } else if (kersel == "uni") {
        Q <- 1.2 
      } 
    } else {
      Q <- 1
    }
  }
  
  # Bias Estimation
  if (g.loo == TRUE) {
    residuals_gs <- ghat_vec.loo - c(ghat)
  } else if (g.loo == FALSE) {
    residuals_gs <- ghat_vec - c(ghat)
  }
  
  Bhat_lp     <- LP.est(y=residuals_gs,x=x,x0=eval,h=h,p=p,v=0,kersel=kersel)
  Bhat_mlp    <- Bhat_lp*Q
  
  # De-biased estimator
  ghat_mlpbc <- ghat - Bhat_mlp
  
  
  #####################  Standard errors and Confidence Intervals  #####################  
  
  # Residuals
  if (res == "loo") {
    if (g.loo == TRUE) {
      epshat <- y - ghat_vec.loo
    } else {
      epshat <- y - ghat_vec
    }
  } else if (res == "cct-hc0" | res == "cct-hc1" | res == "cct-hc2" | res == "cct-hc3") {
    r.pp1 <- matrix(NA, nrow = length(x), ncol = (p + 2))
    for (ip in 1:(p+2)) {
      r.pp1[, ip] <- (x-eval)^(ip-1)
    }
    r.p <- r.pp1[,1:(p+1)]
    K.X    <- (K(u = (x-eval)/h, kersel = kersel)/h)
    invG.q <- qrXXinv((sqrt(K.X)*r.pp1))
    beta.q <- invG.q%*%crossprod(r.pp1*K.X,y)
    
    if (res == "cct-hc0") {
      epshat <- (y - r.pp1%*%beta.q)
    } else {
      Q.vec  <- rowSums((r.pp1 %*% invG.q) * (r.pp1 * K.X))
      if (res == "cct-hc1") {
        epshat <- (y - r.pp1%*%beta.q)/(((length(y) - p+1)/length(y))^(0.5))
      } else if (res == "cct-hc2") {
        epshat <- (y - r.pp1%*%beta.q)/((1-Q.vec)^(0.5))
      } else if (res == "cct-hc3") {
        epshat <- (y - r.pp1%*%beta.q)/((1-Q.vec))
      }
    }
  }
  
  # mPLP Standard errors
  if (fast == 0) {
    se_mplp <- vmPLP(x=x,x0=eval,p=p,h=h,kersel=kersel,epshat = epshat, Q=Q)/sqrt(n*h) 
  } else if (fast == 1) {
    if (kersel == "tri") {
      Kratio <- 0.84
    } else if (kersel == "epa") {
      Kratio <- 0.83 
    } else if (kersel == "uni") {
      Kratio <- 0.86 
    } 
    se_mplp <-  Kratio*(CCT.se(x=x,x0=eval,h=h,p=p,kersel=kersel,epshat = epshat)$se.rb)
  }
  
  
  # Confidence Intervals
  ci_mplp <- c(ghat_mlpbc - qnorm(1-alpha/2)*se_mplp, ghat_mlpbc - qnorm(alpha/2)*se_mplp)
  
  
  # Output
  out <- list(n=n, h=h, Q=Q, g.hat = ghat, g.hat_mlpbc = ghat_mlpbc, 
              se_mplp=se_mplp, ci_mplp=ci_mplp, epshat = epshat)
  
  return(out)
  
}