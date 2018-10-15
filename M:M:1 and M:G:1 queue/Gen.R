# ###############################################################
# #1 compute waiting time
# ###############################################################
# Wait <- function(a,s) #a,s are interarrival & service time
# {
#   m <- length(a)
#   d <- rep(0,m)
#   for (i in 2:m)
#   {
#     d[i] <- max(d[i-1]+s[i-1]-a[i],0)
#   }
#   return(d)
# }
# ############################################

# ################################################################
# #2 generate WT for runlength n, warming up n0
# ################################################################
# WTsample <- function(n,arr,sr,n0)
# {
#   A <- rexp(n+n0,arr)
#   S <- rexp(n+n0,sr)
#   D <- Wait(A,S)[-(1:n0)]
#   return(D)
# }
# ##########################################


###############################################################
#3 AR(1) sample with mu runlength=n, phi1, acvf(0)=r0^2
###############################################################
arsp <- function(n,mu,r0,phi1) #directly simulate from normal to ensure gaussian white noise
{
  Y <- rep(0,n)
  sd1 <- r0*sqrt(1-phi1^2)
  zt <- rnorm(n,0,sd1)
  Y[1] <- rnorm(1,0,r0)
  for (i in 2:n)
  {
    Y[i] <- phi1*Y[i-1] + zt[i]
  }
  return (Y+mu)
}

###############################################################
#4 generating sample path using AR(1) or M/M/1
###############################################################
GenSP <- function(runlen,method,par)
{
  SP <- rep(0,runlen)
  if (method == "mm1")
  {
    # for mm1 par = c( arr, sr, warmingup(=n0))
    # SP <- WTsample(runlen,par[1],par[2],par[3])
    SP <- MM1_sim(par[1],par[2],runlen,par[3])
  }
  else if (method == "ar1")
  {
    # for ar1 par = c( mu, marginal sd(=r0), phi)
    SP <- arsp(runlen,par[1],par[2],par[3])
  }
  else if (method == "mvn")
  {
    # for mvn par = c(mu, marginal sd(=r0), phi)
    SP <- mvrnorm(1,rep(par[1],runlen),par[2]^2*toeplitz(par[3]^(0:(runlen-1))))
  }
  return (SP)
}



