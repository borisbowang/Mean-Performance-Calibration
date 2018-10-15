
################################################################################
# Compute KG given vector a and b
################################################################################

subway <- function(a, b) {  
  #########################################################
  # subway: function to calculate knowledge gradient value
  # 
  # input
  # a: intercept
  # b: slope
  #
  # output:
  # cc: ci values in (24) Qu et al
  # b: bi values in (24) Qu et al
  # idx: nondominated bi indexes
  ###########################################################
  
  K <- length(a)
  ab <- data.frame(a=a, b=b, id=1:K )
  ab <- ab[order(ab$a, decreasing = TRUE), ]
  ab <- ab[!duplicated(ab$b), ]
  K <- nrow(ab)
  ab <- ab[order(ab$b), ]
  cc <- (ab[1,1]-ab[2,1])/(ab[2,2]-ab[1,2])
  idx <- c(1, 2)
  if(K>2) {
    for(i in 2:(K-1)) {
      flag <- 0
      while(flag==0) {
        L <- length(idx)
        M <- L-1
        #compute new intersection 
        cct <- (ab[idx[L],1]-ab[i+1,1])/(ab[i+1,2]-ab[idx[L],2])
        
        #compare the new intersection with the last intersection that has been stored
        if(cct>cc[M]) { #include new intersection
          flag <- 1
          cc <- c(cc,cct)
          idx <- c(idx,i+1)
        } else { #remove the dominated intersection 
          idx <- idx[-L]
          cc <- cc[-M]
          if(length(cc)==0) {
            flag <- 1
            cc <- (ab[idx,1]-ab[i+1,1])/(ab[i+1,2]-ab[idx,2])
            idx <- c(idx,i+1)
          }
        }
        
      }
    }
  }
  return(list(cc=cc, ind=idx, b=ab$b, orig.idx=ab$id))
}



################################################################################
# Given parlist, compute KG as a function of x
################################################################################
Vfun_s <- function(parlist, idx, x)
{
  # idx: selected calibration parameter index
  a_th <- parlist$a[idx]
  c_th <- parlist$c_m[idx, ]
  al_th <- parlist$al[idx]
  be_th <- parlist$be[idx]
  Qi_th <- parlist$Qi_m[idx, , ]

  x_K <- parlist$x_K
  n_K <- length(x_K)
  phi_x <- phi(x, x_K, a_th)

  x_N <- parlist$x_N[!duplicated(parlist$x_N)]
  n_N <- length(x_N)
  xN <- matrix(x_N, n_N, n_K, byrow = F) # N * K
  phi_xN <- apply(xN, 1, function(x){phi(x, x_K, a_th)}) # K * N

  alpha <- t(phi_xN)%*%c_th
  beta <- (t(phi_xN)%*%Qi_th%*%phi_x) * c(sqrt(be_th*(t(phi_x)%*%Qi_th%*%phi_x + 1)/al_th))
  df <- 2*al_th

  # compute
  re <- subway(alpha,beta)
  cc <- re$cc
  ind <- re$ind
  bb <- re$b[ind]
  bdeltavalue <- bb[-1] - bb[-c(length(bb))]
  expectedvalues <- dt(abs(cc), df) * (df + cc^2)/(df - 1)
  expectedvalues <- expectedvalues - abs(cc)*(1 - pt(abs(cc),df))
  kgvalue <- sum(bdeltavalue*expectedvalues)

  return(kgvalue)
}


Vfun_d <- function(parlist, idx, x)
{
  # idx: selected calibration parameter index
  a_th <- parlist$a[idx]
  c_th <- parlist$c_m[idx, ]
  al_th <- parlist$al[idx]
  be_th <- parlist$be[idx]
  Qi_th <- parlist$Qi_m[idx, , ]
  Ri_th <- parlist$Ri_m[idx, , ]
  z1 <- sum(Ri_th)

  x_K <- parlist$x_K
  n_K <- length(x_K)
  phi_x <- phi(x, x_K, a_th)

  x_N <- parlist$x_N[!duplicated(parlist$x_N)]
  n_N <- length(x_N)
  xN <- matrix(x_N, n_N, n_K, byrow = F) # N * K
  phi_xN <- apply(xN, 1, function(x){phi(x, x_K, a_th)}) # K * N

  alpha <- t(phi_xN)%*%c_th
  beta <- (t(phi_xN)%*%Qi_th%*%phi_x) * c(sqrt(be_th*z1*(t(phi_x)%*%Qi_th%*%phi_x*z1 + 1)/al_th))
  df <- 2*al_th

  # compute
  re <- subway(alpha,beta)
  cc <- re$cc
  ind <- re$ind
  bb <- re$b[ind]
  bdeltavalue <- bb[-1] - bb[-c(length(bb))]
  expectedvalues <- dt(abs(cc), df) * (df + cc^2)/(df - 1)
  expectedvalues <- expectedvalues - abs(cc)*(1 - pt(abs(cc),df))
  kgvalue <- sum(bdeltavalue*expectedvalues)

  return(kgvalue)
}



################################################################################
# Select decision max kg value, optimize over x
################################################################################
dec_s <- function(parlist, idx, lob_x, upb_x)
{
  par <- 1
  kgfun <- function(par) -Vfun_s(parlist, idx, par)
  res <- optim(par, kgfun, method = "Brent", lower = lob_x, upper = upb_x)
  return(res$par)
}


dec_d <- function(parlist, idx, lob_x, upb_x)
{
  par <- 1
  kgfun <- function(par) -Vfun_d(parlist, idx, par)
  res <- optim(par, kgfun, method = "Brent", lower = lob_x, upper = upb_x)
  return(res$par)
}



################################################################################
# Select optimal calibration parameter theta & decision x
################################################################################
fdec <- function(parlist, target, lob_x, upb_x)
{
  # target matrix = n_th * n_K (with same rows)
  mse <- apply((parlist$post_mu-target)^2 + parlist$post_va, 1, sum)
  idx <- which.min(mse)

  par <- (lob_x + upb_x)/2
  x_K <- parlist$x_K
  a_th <- parlist$a[idx]
  c_th <- parlist$c_m[idx, ]

  objfun <- function(par) -t(phi(par, x_K, a_th))%*%c_th
  res <- optim(par, objfun, method = "Brent", lower = lob_x, upper = upb_x)
  return(list(x = res$par, th = idx))
}






