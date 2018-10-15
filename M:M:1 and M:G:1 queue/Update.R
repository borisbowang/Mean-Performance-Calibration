#################################################################
#1.1 functions of Local Time Method (LTM) T&S
#################################################################
#(a) compute local time
ltime <- function(theta,sig2,tar=0) #tar is target c
{
  z <- -abs(theta-tar)/sqrt(sig2)
  lt <- 2*sqrt(sig2)*(z*pnorm(z)+dnorm(z))
  return(lt)
}

#################################################################
#1.2 rho estimation
#################################################################
rho_hat <- function(Y)
{
  X <- Y - mean(Y)
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
#2 One Step Update
#################################################################
updates <- function(parlist,idx,dat)
{
  #parlist = list(theta_s,tau,a,b,var_s)
  #dat = ybar
  a <- parlist$a[idx]
  b <- parlist$b[idx]
  mu <- parlist$theta_s[idx]
  tau <- parlist$tau[idx]

  parlist$a[idx] <- a + 1/2
  parlist$theta_s[idx] <- mu +(dat-mu)/(tau + 1)
  parlist$b[idx] <- b + tau*(dat-mu)^2/(2*(tau + 1))
  parlist$tau[idx] <- tau + 1
  parlist$var_s[idx] <- parlist$b[idx]/(parlist$tau[idx]*(parlist$a[idx]-1))

  return (parlist)
}


updated_kn <- function(parlist,idx,Rr,dat)
{
  #parlist = list(theta_d,q,al,be,var_d)
  #dat = Y
  L <- length(dat)
  v <- rep(1,L)
  mu <- parlist$theta_d[idx]
  q <- parlist$q[idx]
  al <- parlist$al[idx]
  be <- parlist$be[idx]

  parlist$al[idx] <- al + L/2
  parlist$theta_d[idx] <- mu + (t(v)%*%Rr%*%(dat-mu))/(q+sum(Rr))
  # parlist$be[idx] <- be + (q*t(dat-mu)%*%Rr%*%(dat-mu))/(2*(q+sum(Rr)))
  parlist$be[idx] <- be + (q*t(dat-mu)%*%Rr%*%(dat-mu)+t(v)%*%Rr%*%(v%*%t(dat)-dat%*%t(v))%*%Rr%*%dat)/(2*(q+sum(Rr)))
  # parlist$be[idx] <- be + (t(dat)%*%Rr%*%dat)/2 + q*mu^2/2 - (q*mu+t(v)%*%Rr%*%dat)^2/(2*(q+sum(Rr)))
  parlist$q[idx] <- q+sum(Rr)
  parlist$var_d[idx] <- parlist$be[idx]/(parlist$q[idx]*(parlist$al[idx]-1))
  
  return (parlist)
}


updated_sw <- function(parlist,idx,dat)
{
  #parlist = list(theta_d,q,Ydat,al,be,var_d)
  #Ydat is also a list
  #dat = Y
  L <- length(dat)
  v <- rep(1,L)
  Ypaths <- rbind(parlist$Ydat[[idx]],dat)
  Rr <- RR(rho_hat(Ypaths),L)
  mu <- parlist$theta_d[idx]
  q <- parlist$q[idx]
  al <- parlist$al[idx]
  be <- parlist$be[idx]

  parlist$al[idx] <- al + L/2
  parlist$theta_d[idx] <- mu + (t(v)%*%Rr%*%(dat-mu))/(q+sum(Rr))
  # parlist$be[idx] <- be + (q*t(dat-mu)%*%Rr%*%(dat-mu))/(2*(q+sum(Rr)))
  parlist$be[idx] <- be + (q*t(dat-mu)%*%Rr%*%(dat-mu)+t(v)%*%Rr%*%(v%*%t(dat)-dat%*%t(v))%*%Rr%*%dat)/(2*(q+sum(Rr)))
  # parlist$be[idx] <- be + (t(dat)%*%Rr%*%dat)/2 + q*mu^2/2 - (q*mu+t(v)%*%Rr%*%dat)^2/(2*(q+sum(Rr)))
  parlist$q[idx] <- q+sum(Rr)
  parlist$var_d[idx] <- parlist$be[idx]/(parlist$q[idx]*(parlist$al[idx]-1))
  parlist$Ydat[[idx]] <- Ypaths
  
  return (parlist)
}









