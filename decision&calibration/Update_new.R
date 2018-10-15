#################################################################
#1.0 phi vector: set basis function (added 07/24/18)
#################################################################
phi <- function(x, x_K, a = 0.5)
{
  K <- length(x_K)
  phi_x <- exp(-(x - x_K)^2*a)
  
  return (phi_x)
}


#################################################################
#1.1 functions of Local Time Method (LTM) T&S
#################################################################
#(a) compute local time
ltime <- function(theta,sig2,tar=0) 
{
  # tar is target c
  # theta, sig2 are posterior mean and variance
  z <- -abs(theta-tar)/sqrt(sig2)
  lt <- 2*sqrt(sig2)*(z*pnorm(z)+dnorm(z))
  return (lt)
}

#################################################################
#1.2 rho estimation (modified 09/08/18)
#################################################################
rho_hat <- function(Y)
{
  # X <- Y - mean(Y)
  # Y is collection of de-meaned sample paths
  X <- Y
  if ( is.null(nrow(X)) || nrow(X)==0 ){
    L <- length(X)
    r_h <- sum(X[-1]*X[-L])/sum(X[2:(L-1)]^2)
  }
  else{
    L <- ncol(X)
    r_h <- sum(X[,-1]*X[,-L])/sum(X[,2:(L-1)]^2)
  }
  return (min(r_h*(L-2)/(L-1),0.999999))
}

#################################################################
#1.3 rho -> Rr (reversed corr matrix)
#################################################################
RR <- function(rho,runlen)
{
  Rr <- toeplitz(c(1+rho^2 , -rho, rep(0,runlen-2)))
  Rr[1,1] <- Rr[runlen,runlen] <- 1
  return ( Rr/(1-rho^2) )
}

#################################################################
#1.4 random select index (when there is tie)
#################################################################
r_select <- function(x){
  if (length(x) == 1){
    return (x)
  }else{
    return (sample(x,1))
  }
}

  

#################################################################
#2 One Step Update with matrix para settings 
#################################################################
update_s <- function(parlist, idx, x, ybar)
{
  # update for summary stat
  # parlist = list(x_K, a, c_m, Q_m, Qi_m, al, be, post_mu, post_va, x_N) # m: method
  # idx = calibration parameter index selected
  # x = decision selected
  parlist$x_N <- c(parlist$x_N, x)

  x_K <- parlist$x_K
  a_th <- parlist$a[idx] # a related to theta
  phi_x <- phi(x, x_K, a_th)

  c_th <- parlist$c_m[idx, ]
  Q_th <- parlist$Q_m[idx, , ]
  al_th <- parlist$al[idx]
  be_th <- parlist$be[idx]
  Qi_th <- parlist$Qi_m[idx, , ] # inverse

  parlist$Q_m[idx, , ] <- Q_new <- Q_th + phi_x%*%t(phi_x)
  parlist$Qi_m[idx, , ] <- Qi_new <- solve(Q_new)
  parlist$c_m[idx, ] <- c_new <- c( Qi_new%*%(Q_th%*%c_th + ybar*phi_x) )
  parlist$al[idx] <- al_new <- al_th + 1/2
  parlist$be[idx] <- be_new <- be_th + (ybar^2 + t(c_th)%*%Q_th%*%c_th - t(c_new)%*%Q_new%*%c_new)/2
  

  n_K <- length(x_K)
  xK <- matrix(x_K, n_K, n_K, byrow = F)
  phi_xK <- apply(xK, 1, function(x){return(phi(x, x_K, a_th))}) # symmetric (=Sigma K*K)

  parlist$post_mu[idx, ] <- phi_xK%*%c_new
  parlist$post_va[idx, ] <- diag(phi_xK%*%Qi_new%*%phi_xK) * c(be_new/(al_new - 1))
  
  return(parlist)
}


# update_d <- function(parlist, idx, x, Y)
# {
#   # update for detailed method
#   # parlist = list(x_K, a, Ydat, c_m, Q_m, Qi_m, al, be, post_mu, post_va, x_N, Ri_m)
#   # Ydat is list, each Y[[idx]] is matrix containing de-meaned sample paths (rbind)
#   # Y is NOT de-meaned
#   parlist$x_N <- c(parlist$x_N, x)

#   x_K <- parlist$x_K
#   a_th <- parlist$a[idx]
#   phi_x <- phi(x, x_K, a_th)

#   L <- length(Y)
#   v1 <- rep(1, L) # eqiv as as.matrix(rep(1,L))
#   Ys <- rbind(parlist$Ydat[[idx]], Y - mean(Y)) # de-meaned sample path
#   Ri <- RR(rho_hat(Ys),L) # inverse i.e. Omega in paper
#   z1 <- t(v1)%*%Ri%*%v1 # equiv sum(Ri)
#   z2 <- t(v1)%*%Ri%*%Y
#   z3 <- t(Y)%*%Ri%*%Y

#   c_th <- parlist$c_m[idx, ]
#   Q_th <- parlist$Q_m[idx, , ]
#   al_th <- parlist$al[idx]
#   be_th <- parlist$be[idx]
#   Qi_th <- parlist$Qi_m[idx, , ] # inverse

#   parlist$Q_m[idx, , ] <- Q_new <- Q_th + c(z1)*( phi_x%*%t(phi_x) )
#   parlist$Qi_m[idx, , ] <- Qi_new <- solve(Q_new)
#   parlist$c_m[idx, ] <- c_new <- c( Qi_new%*%(Q_th%*%c_th + c(z2)*phi_x) )
#   parlist$al[idx] <- al_new <- al_th + L/2
#   parlist$be[idx] <- be_new <- be_th + (z3 + t(c_th)%*%Q_th%*%c_th - t(c_new)%*%Q_new%*%c_new)/2
  

#   parlist$Ydat[[idx]] <- Ys
#   parlist$Ri_m[idx, , ] <- Ri

#   n_K <- length(x_K)
#   xK <- matrix(x_K, n_K, n_K, byrow = F)
#   phi_xK <- apply(xK, 1, function(x){return(phi(x, x_K, a_th))}) # symmetric (=Sigma K*K)

#   parlist$post_mu[idx, ] <- phi_xK%*%c_new
#   parlist$post_va[idx, ] <- diag(phi_xK%*%Qi_new%*%phi_xK) * c(be_new/(al_new - 1))

#   return (parlist)
# }

update_d <- function(parlist, idx, x, Y)
{
  # update for detailed method (with known Ri_m)
  # parlist = list(x_K, a, c_m, Q_m, Qi_m, al, be, post_mu, post_va, x_N, Ri_m)
  
  parlist$x_N <- c(parlist$x_N, x)

  x_K <- parlist$x_K
  a_th <- parlist$a[idx]
  phi_x <- phi(x, x_K, a_th)

  L <- length(Y)
  v1 <- rep(1, L) # eqiv as as.matrix(rep(1,L))

  Ri <- parlist$Ri_m[idx, , ] # inverse i.e. Omega in paper
  z1 <- t(v1)%*%Ri%*%v1 # equiv sum(Ri)
  z2 <- t(v1)%*%Ri%*%Y
  z3 <- t(Y)%*%Ri%*%Y

  c_th <- parlist$c_m[idx, ]
  Q_th <- parlist$Q_m[idx, , ]
  al_th <- parlist$al[idx]
  be_th <- parlist$be[idx]
  Qi_th <- parlist$Qi_m[idx, , ] # inverse

  parlist$Q_m[idx, , ] <- Q_new <- Q_th + c(z1)*( phi_x%*%t(phi_x) )
  parlist$Qi_m[idx, , ] <- Qi_new <- solve(Q_new)
  parlist$c_m[idx, ] <- c_new <- c( Qi_new%*%(Q_th%*%c_th + c(z2)*phi_x) )
  parlist$al[idx] <- al_new <- al_th + L/2
  parlist$be[idx] <- be_new <- be_th + (z3 + t(c_th)%*%Q_th%*%c_th - t(c_new)%*%Q_new%*%c_new)/2
  

  n_K <- length(x_K)
  xK <- matrix(x_K, n_K, n_K, byrow = F)
  phi_xK <- apply(xK, 1, function(x){return(phi(x, x_K, a_th))}) # symmetric (=Sigma K*K)

  parlist$post_mu[idx, ] <- phi_xK%*%c_new
  parlist$post_va[idx, ] <- diag(phi_xK%*%Qi_new%*%phi_xK) * c(be_new/(al_new - 1))

  return(parlist)
}



#######################################################################
#3 Compute posterior mean & var of expected output at x_K from parlist
#######################################################################
postmv <- function(parlist)
{
  # initially compute
  c_m <- parlist$c_m
  n_th <- nrow(c_m)
  n_K <- ncol(c_m)
  post_mu <- post_va <- matrix(0, n_th, n_K)

  x_K <- parlist$x_K
  a_m <- parlist$a
  al <- parlist$al
  be <- parlist$be
  Qi_m <- parlist$Qi_m # inverse

  for (i in 1:n_th){
    a_th <- a_m[i]
    c_th <- c_m[i, ]
    Qi_th <- Qi_m[i, , ]
    al_th <- al[i]
    be_th <- be[i]

    xK <- matrix(x_K, n_K, n_K, byrow = F)
    phi_xK <- apply(xK, 1, function(x){return(phi(x, x_K, a_th))}) # symmetric (=Sigma K*K)

    post_mu[i, ] <- phi_xK%*%c_th
    post_va[i, ] <- diag(phi_xK%*%Qi_th%*%phi_xK) * c(be_th/(al_th - 1))
  }

  return(list(post_mu = post_mu, post_va = post_va))
}











