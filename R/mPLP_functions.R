######################################################################################
## Auxiliary functions for the R programs mPLPnp and mPLPrd
## Authors: Giuseppe Cavaliere, Sílvia Gonçalves, Morten Ø. Nielsen, & Edoardo Zanelli
## version: 1.01 (9.12.2025)
######################################################################################

library(sandwich)

# Kernel function 
K <- function(u, kersel) {
  
  result <- if (is.character(kersel) && length(kersel) == 1 &&
                kersel %in% c("epanechnikov", "epa")) {
    
    ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0)
    
  } else if (is.character(kersel) && length(kersel) == 1 &&
             kersel %in% c("triangular", "tri")) {
    
    ifelse(abs(u) <= 1, (1 - abs(u)), 0)
    
  } else if (is.character(kersel) && length(kersel) == 1 &&
             kersel %in% c("uniform", "uni")) {
    
    ifelse(abs(u) <= 1, 0.5, 0)
    
  } else {
    stop("Invalid kernel specification.")
  }
  
  return(as.numeric(result))  
}


# Polynomial basis function
r_p <- function(p) {
  function(u) {
    u <- as.vector(u)
    outer(u, 0:p, `^`)
  }
}


# Inverse of the Gamma_p matrix
Gamma.p.inv <- function(u,h,p,kersel) {
  n <- length(u)
  X <- matrix(NA, nrow = length(u), ncol = p+1)
  for (i in 1:(p+1)) {
    X[,i] <- (u)^(i-1)
  }
  K.X <- K(u = u, kersel = kersel)/h
  nonzero_indices <- which(K.X != 0) 
  X <- X[nonzero_indices, ]
  K.X <- K.X[nonzero_indices]
  Gamma.inv <- qrXXinv((sqrt(K.X)*X))*n
}


# Local polynomial estimator 
LP.est <- function(y,x,x0,h,p,v,kersel) {
  
  n <- length(x)
  X <- matrix(NA, nrow = n, ncol = p+1)
  u <- (x-x0)/h
  for (i in 1:(p+1)) {
    X[,i] <- (u)^(i-1)
  }
  K.X <- K(u = u, kersel = kersel)/h
  nonzero_indices <- which(K.X != 0) 
  X <- X[nonzero_indices, ]
  K.X <- K.X[nonzero_indices]
  
  bwcheck <- 21
  if (!is.null(bwcheck)) {
    bw.min   <- sort(abs(x-x0))[bwcheck]
    h     <- max(h, bw.min)
  }
  
  Gamma.inv <- Gamma.p.inv(u=u, h=h, p=p, kersel=kersel) 
  
  
  e1 <- numeric(p+1)
  e1[v+1] <- 1
  ghat <- factorial(v)*e1%*%Gamma.inv%*%crossprod(X*K.X,y[nonzero_indices])/(n*(h^v))
  return(ghat)
}


# Equivalent kernel of a local polynomial estimator evaluated at a single x_i  
wix0 <- function(x,xi,x0,h,p,v,kersel) {
  bwcheck <- 21
  if (!is.null(bwcheck)) {
    bw.min   <- sort(abs(x-x0))[bwcheck]
    h        <- max(h, bw.min)
  }
  e1 <- numeric(p+1)
  e1[v+1] <- 1
  n <- length(x)
  u <- (x-x0)/h
  Gamma.inv <- Gamma.p.inv(u=u, h=h, p=p, kersel=kersel)
  
  W <- numeric(length(xi))
  u_all <- (xi - x0) / h
  K_all <- K(u = u_all, kersel = kersel)
  nonzero_indices.2 <- which(K_all != 0)
  for (ii in nonzero_indices.2) {
    ui <- u_all[ii]
    kval <- K_all[ii]
    ru <- kval * ui^(0:p)
    W[ii] <- factorial(v) * e1 %*% Gamma.inv %*% ru
  }
  return(W)
}

## -- Equivalent kernel of a (mPLP based) debiased LP estimator evaluated at a single x_i -- ## 

wmPLP <- function(x,xi,x0,h,p,kersel, Q) {
  n <- length(x)
  A <- c(1+Q)*wix0(x = x,xi = xi, x0 = x0, h = h, p=p, v=0, kersel = kersel) 
  B <- numeric(n)
  
  U.Z = abs((x-x0)/h) <= 1
  nonzero_indices <- which(U.Z != 0) 
  for (jj in nonzero_indices) {
    BB <- c(Q) * wix0(x=x, xi=x[jj], x0=x0, h=h, p=p, v=0, kersel=kersel)*wix0(x=x, xi=xi, x0=x[jj], h=h, p=p, v=0, kersel=kersel)/(n*h)
    B[jj] <- BB
  } 
  B <- sum(B)
  A - B
}


# Variance estimator of a (mPLP-based) bias corrected test statistic 
vmPLP <- function(x,x0,h,p,kersel,epshat,Q) {
  n <- length(x)
  C <- numeric(n)
  U.Z <- (abs((x-x0)/(2*h)) <= 1) 
  nonzero_indices <- which(U.Z != 0) 
  for (ix in seq_along(nonzero_indices)) {
    ii <- nonzero_indices[ix]
    C[ii] <- ((wmPLP(x = x, xi = x[ii], x0 = x0, h=h, p=p, kersel=kersel, Q=Q)*epshat[ii])^2)/(n*h)
  } 
  C <- sum(C)
  sqrt(C)
}


# CCT (single bandwidth) standard errors  
CCT.se <- function(x,x0,h,p,kersel,epshat) {
  
  n <- length(x)
  r.pp1 <- matrix(NA, nrow = n, ncol = (p + 2))
  for (ip in 1:(p+2)) {
    r.pp1[, ip] <- (x-x0)^(ip-1)
  }
  r.p <- r.pp1[,1:(p+1)]
  
  e.p1    <- matrix(0,(p+2),1); e.p1[p+2]=1
  K.X <- K(u = (x-x0)/h, kersel = kersel)/h
  invG.q  <- qrXXinv((sqrt(K.X)*r.pp1))
  invG.p  <- qrXXinv((sqrt(K.X)*r.p))
  L <- crossprod(r.p*K.X,((x-x0)/h)^(p+1))
  
  Q.q <- t(t(r.p*K.X) - (h^(p+1))*(L%*%t(e.p1))%*%t(t(invG.q%*%t(r.pp1))*K.X))
  Q.p <- r.p*K.X
  
  MM.p <- crossprod(c(epshat)*Q.p)
  MM.q <- crossprod(c(epshat)*Q.q)
  MMM.p <- invG.p%*%MM.p%*%invG.p
  MMM.q <- invG.p%*%MM.q%*%invG.p
  
  
  M.p <- sqrt(MMM.p[1,1])
  M.q <- sqrt(MMM.q[1,1])
  
  ses <- list(se.us = M.p, se.rb=M.q)
  return(ses)
}

qrXXinv = function(x, ...) {
  chol2inv(chol(crossprod(x)))
}