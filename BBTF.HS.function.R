library(LaplacesDemon)
library(mvtnorm)
library(GIGrvg)
library(coda)
library(pgdraw)
library(tmg)



## ------------------------------------- ##
#              BBTF horseshoe             #
## ------------------------------------- ##


BBTF.HS <- function(y, x=1:length(y), mc=10000, burn=2000, eta=500, k=1, upper=T, init.sig=c(0.1,0.1), init.u=c(1,1), shape="N", itr=5000){
  
  if(shape != "N" & shape != "NI" & shape != "ND"){stop("shape: N (none), NI (nearly isotonic) or ND (nearly decreasing) ")}
  
  trunc <- 10^(-6)
  n <- length(y)
  
  # hyper parameter
  a_sig <- init.sig[1]
  b_sig <- init.sig[2]
  a_u <- init.u[1]
  b_u <- init.u[2]
  
  # difference operator
  D1 <- matrix(0, n-1,n)
  for (i in 1:n-1) {
    for (j in 1:n) {
      if(i == j){
        D1[i,j] <- 1
      }
      if(i+1 == j){
        D1[i,j] <- -1
      }
    }
  }
  if(k == 0){
    DD <- D1
    D <- rbind(c(1,rep(0, n-1)),DD)
  }
  if(k == 1){
    delta1 <- diff(x)
    DD <- D1[1:(n-2),1:(n-1)]%*%diag(k/delta1)%*%D1
    D <- rbind(t(matrix(c(1,rep(0, n-1), 0,1,rep(0, n-2)),n, 2)),DD)
  }
  if(k == 2){
    DD1 <- D1[1:(n-2),1:(n-1)]%*%diag(1/diff(x))%*%D1
    DD2 <- D1[1:(n-3),1:(n-2)]%*%diag(2/diff(x,lag=2))%*%DD1
    D <- rbind(t(matrix(c(1,rep(0, n-1), 0,1,rep(0, n-2), 0,0, 1, rep(0, n-3)),n, 3)),DD2)
  }
  
  TD <- t(D)
  TD1 <- t(D1)
  
  # parameter
  theta.pos <- matrix(0, mc, n)
  sig.pos <- numeric(mc)
  
  if(shape == "N"){
    
    ## no constraint
    print("shape: none")
    
    # hyper parameter
    kappa <- rep(1/2, n)
    
    # initial value
    xi <- rep(1,n)
    omega <- rep(1,n)
    sig <- 1
    u <- rep(1, n)
    tau <- 1
    tauvec <- c(rep(1, 1+k), rep(tau, n-k-1))
    nu <- rep(1, n-k-1)
    psi <- 1
    
    for (i in 1:mc) {
      
      # sample omega
      omega <- pgdraw(rep(1,n), eta*xi)
      
      # sample theta
      TDUD <- TD%*%diag(1/(tauvec*u))%*%D
      A <- diag(1/sig,n) + TDUD/sig
      if(upper==FALSE){
        b <- TDUD%*%y/sig
      }else{
        b <- -TDUD%*%y/sig
      }
      Omega <- diag(omega,n)
      Sigma <- solve(eta^2*Omega + A)
      Sigma <- (t(Sigma)+Sigma)/2
      mu <- Sigma%*%(eta*kappa + b)
      xi <- as.vector(rmvnorm(1, mu, Sigma))
      
      if(upper==FALSE){
        th <- y-xi
      }else{
        th <- xi+y
      }
      theta.pos[i,] <- th
      Dth <- D%*%th
      
      # sample sigma^2
      aa <- t(th)%*%TDUD%*%th/2 + t(xi)%*%(xi)/2 + b_sig
      sig <- 1/rgamma(1, a_sig+n, aa )
      sig.pos[i] <- sig
      
      # sample u^2
      aa <- (D%*%th)^2/(2*sig*tauvec)
      u[1:(k+1)] <- 1/rgamma(k+1, 1/2+a_u, aa[1:(k+1)]+b_u)
      u[(k+2):n] <- 1/rgamma(n-k-1, 1, aa[(k+2):n]+1/nu)
      u[u<trunc] <- trunc
      
      # sample tau^2
      aa <-  sum( ((D%*%th)^2/(2*sig*u))[(k+2):n] ) + 1/psi
      tau <- 1/rgamma(1, (n-k)/2, aa)
      tau[tau<trunc] <- trunc
      tauvec <- c(rep(1, 1+k), rep(tau, n-k-1))
      
      # sample nu
      nu <- 1/rgamma(n-k-1, 1/2, 1/u[(k+2):n]+1)
      
      # sample psi
      psi <- 1/rgamma(1, 1/2, 1/tau+1)
      
      if(!is.na(itr)){
        if(i/ itr == round(i/itr)){print(i)}
      }
    }
  }else if(shape == "NI" | shape == "ND"){
    
    if(shape == "NI"){
      print("shape: nearly isotonic")
    }
    if(shape == "ND"){
      print("shape: nearly decreasing")
      D1 <- -D1
      TD1 <- t(D1)
    }
    
    # hyper parameter
    kappa <- rep(1/2, n)
    init.rho=c(1,1)
    a_rho <- init.rho[1]
    b_rho <- init.rho[2]
    
    # initial value
    xi <- rep(1,n)
    omega <- rep(1,n)
    sig <- 1
    u <- rep(1, n)
    tau <- 1
    tauvec <- c(rep(1, 1+k), rep(tau, n-k-1))
    nu <- rep(1, n-k-1)
    psi <- 1
    vv <- rep(1, n-1)
    rho <- 1
    
    for (i in 1:mc) {
      
      # sample omega
      omega <- pgdraw(rep(1,n), eta*xi)
      
      # sample theta
      TDUD <- TD%*%diag(1/(tauvec*u))%*%D
      TD1VD1 <- TD1%*%diag(1/(2*vv*rho))%*%D1/sig
      A <- diag(1/sig,n) + TDUD/sig + TD1VD1 
      if(upper==FALSE){
        b <- TDUD%*%y/sig + TD1VD1%*%y +TD1%*%rep(1/(2*sig*rho), n-1)
      }else{
        b <- -(TDUD%*%y/sig + TD1VD1%*%y + TD1%*%rep(1/(2*sig*rho), n-1) )
      }
      Omega <- diag(omega,n)
      Sigma <- solve(eta^2*Omega + A)
      Sigma <- (t(Sigma)+Sigma)/2
      mu <- Sigma%*%(eta*kappa + b)
      xi <- as.vector(rmvnorm(1, mu, Sigma))
      
      if(upper==FALSE){
        th <- y-xi
      }else{
        th <- xi+y
      }
      theta.pos[i,] <- th
      D1th <- D1%*%th
      
      # sample sigma^2
      aa <- t(th)%*%TDUD%*%th/2 + t(xi)%*%(xi)/2 + b_sig + sum(D1th^2/(4*vv*rho))
      sig <- 1/rgamma(1, a_sig+(3*n-1)/2, aa )
      sig.pos[i] <- sig
      
      # sample u^2
      aa <- (D%*%th)^2/(2*sig*tauvec)
      u[1:(k+1)] <- 1/rgamma(k+1, 1/2+a_u, aa[1:(k+1)]+b_u)
      u[(k+2):n] <- 1/rgamma(n-k-1, 1, aa[(k+2):n]+1/nu)
      u[u<trunc] <- trunc
      
      # sample tau^2
      aa <-  sum( ((D%*%th)^2/(2*sig*u))[(k+2):n] ) + 1/psi
      tau <- 1/rgamma(1, (n-k)/2, aa)
      tau[tau<trunc] <- trunc
      tauvec <- c(rep(1, 1+k), rep(tau, n-k-1))
      
      # sample nu
      nu <- 1/rgamma(n-k-1, 1/2, 1/u[(k+2):n]+1)
      
      # sample psi
      psi <- 1/rgamma(1, 1/2, 1/tau+1)
      
      # sample v
      vv <- as.vector(tapply(D1th^2/(2*rho*sig), 1:(n-1), rgig, n=1, lambda=1/2, psi=1/(2*rho*sig)))
      
      # sample rho
      rho <- 1/rgamma(1, (n-1)/2+a_rho, sum(D1th^2/(4*vv*sig)) + b_rho)
      
      if(!is.na(itr)){
        if(i/ itr == round(i/itr)){print(i)}
      }
    }
  }
  
  om <- 1:burn
  result <- list(theta.pos=theta.pos[-om,], sig.pos=sig.pos[-om])
  return(result)
  
}











