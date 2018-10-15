#####################################
# author: Bo Wang
# latest update: 09/17/2018
#####################################
rm(list=ls())
source("QN_sim.R")
source("Update_new.R")
source("KGfun.R")

require(mlegp)
# require(ggplot2)
# require(reshape2)

set.seed(11) #set seed

###################################################################################
# 1 Initialization and parameters settings
###################################################################################
# runlength
runlen <- 50
# warmup (1000)
warmup <- 1000

# simulation budget
N <- 1000
# micro replications (100)
IT <- 100


# arrival rate = 0.5
parA <- 0.5
# True probability
pr <- rbind(c(0, 0.8, 0.2, 0),
            c(0.2, 0, 0, 0.8),
            c(0.2, 0, 0, 0.8))


parS3 <- 1
# decision variable 
# service rate 1,2
lob_x <- 0.8
upb_x <- 1.4


m_cyc <- function(parS1, p) # compute true mean cycle time
{
  arr1 <- parA / (1 - p[1, 2] * p[2, 1] - p[1, 3]*p[3, 1])
  arr2 <- arr1 * p[1, 2]
  arr3 <- arr1 * p[1, 3]

  parS2 <- 2 - parS1
  EN1 <- arr1 / (parS1 - arr1)
  EN2 <- arr2 / (parS2 - arr2)
  EN3 <- arr3 / (parS3 - arr3)
  cyctime <- (EN1 + EN2 + EN3) / parA
  return(cyctime)
}


# compute true optimal decision ########################
x_opt <- optim(par=1, function(x){m_cyc(x, pr)}, method = "Brent", 
               lower = lob_x, upper = upb_x)$par
mt_opt <- m_cyc(x_opt, pr)

# observed decision points #############################
# x_K <- seq(0.8, 1.4, by=0.2)
x_K <- seq(0.8, 1.4, length.out = 7) # adjust K
n_K <- length(x_K)
mt_K <- m_cyc(x_K, pr)


# calibration alternatives #############################
# p12 <- seq(0.7, 0.9, by=0.02)
p12 <- seq(0.75, 0.85, length.out = 11) # adjust # of calibration alternatives
# p12 <- c(0.6, 0.8)  # 2 calibration alternatives
n_th <- length(p12)
target <- matrix(-mt_K, n_th, n_K, byrow = T) # target matrix
pr_m <- array(0, dim=c(n_th, 3, 4))
mt_m <- matrix(0, n_th, n_K)

for (i in 1:n_th){
  pr_m[i, , ] <- rbind(c(0, p12[i], (1 - p12[i]), 0),
                        c(0.2, 0, 0, 0.8),
                        c(0.2, 0, 0, 0.8))
  mt_m[i, ] <- -m_cyc(x_K, pr_m[i, , ])
}
true.mse <- apply((mt_m - target)^2, 1, sum)
trueid <- which.min(true.mse)



###################################################################################
# 2 Initial parameters, priors
###################################################################################

# direct setting prior #############################
al <- rep(3/2, n_th)
be <- runif(n_th, 200, 500)

post_mu <- post_va <- matrix(0, n_th, n_K)


Ri_m <- array(0, dim=c(n_th, runlen, runlen))
for (i in 1:n_th){
  rho <- rep(0, n_K)
  for (j in 1:n_K){
    yt <- as.ts(-QN_sim(parA, c(x_K[j], (2-x_K[j]), parS3), pr_m[i, , ], 5000, warmup))
    rho[j] <- c(acf(yt, lag.max=1, type="correlation", plot=F)$acf)[2]
  }
  rho_avg <- mean(rho)
  Ri_m[i, , ] <- RR(rho_avg, runlen)
}


# para/prior for fix rank kriging ##################
# use some initial runs ############################

a <- rep(0.5, n_th)
c_m <- matrix(rnorm(n_th*n_K), n_th, n_K)
Q_m <- Qi_m <- array(0, dim=c(n_th, n_K, n_K))

xK <- matrix(x_K, n_K, n_K, byrow = F)


# use true mean response ##########
# y_gp <- mlegp(x_K, -mt_K, constantMean=1, nugget.known=0)
# rr <- y_gp$beta

for (i in 1:n_th){

  y_gp <- mlegp(x_K, -mt_m[i, ], constantMean=1, nugget.known=0)
  a[i] <- y_gp$beta
  # y_gp <- km(design = as.matrix(x_K), response = -mt_m[i, ], covtype="gauss")
  # a[i] <- 1/(y_gp@covariance@range.val^2*2)
  while (a[i] > 200){
    y_gp <- mlegp(x_K, -mt_m[i, ], constantMean=1, nugget.known=0)
    a[i] <- y_gp$beta
  }
  # a[i] <- rr
  Qi_m[i, , ] <- phi_xK <- apply(xK, 1, function(x){return(phi(x, x_K, a[i]))})
  Q_m[i, , ] <- solve(phi_xK)
  # Qi_m[i, , ] <- toeplitz(0^(0:(n_K-1)))
  # Q_m[i, , ] <- RR(0, n_K)
  c_m[i, ] <- Q_m[i, , ]%*%mt_m[i, ]

}


for (it in 1:IT){

  # n_0 <- 10
  # for (i in 1:n_th){
  #   y <- matrix(0, n_0, n_K)
  #   for (j in 1:n_K){
  #     for (k in 1:n_0){
  #       y[k, j] <- mean(-QN_sim(parA, c(x_K[j], (2-x_K[j]), parS3), pr_m[i, , ], runlen, warmup))
  #     }
  #   }
  #   ybars <- apply(y, 2, mean)
  #   yvars <- apply(y, 2, var)/n_0
  #   y_gp <- mlegp(x_K, ybars, constantMean=1, nugget=yvars, nugget.known=1)

  #   a[i] <- y_gp$beta
  #   # Qi_m[i, , ] <- phi_xK <- apply(xK, 1, function(x){return(phi(x, x_K, a[i]))}) + yvars*diag(rep(1,n_K))/(y_gp$sig2)
  #   # Q_m[i, , ] <- solve(phi_xK)
  #   Qi_m[i, , ] <- phi_xK <- apply(xK, 1, function(x){return(phi(x, x_K, a[i]))})
  #   Q_m[i, , ] <- solve(phi_xK)
  #   c_m[i, ] <- Q_m[i, , ]%*%ybars

  # }


  # parameter list for 2 methods #####################
  parlist_s <- list(x_K = x_K, a = a, c_m = c_m, Q_m = Q_m, Qi_m = Qi_m, al = al, 
                    be = be, post_mu = post_mu, post_va = post_va, x_N = x_K)

  parlist_d <- list(x_K = x_K, a = a, c_m = c_m, Q_m = Q_m, Qi_m = Qi_m, al = al, 
                    be = be, post_mu = post_mu, post_va = post_va, x_N = x_K, Ri_m = Ri_m)



  ###################################################################################
  # 3 Initial runs
  ###################################################################################

  # num of initial runs for each calibration alternative (>= K) ########
  n_init <- 3
  # x_init <- c(x_K, runif((n_init - n_K), lob_x, upb_x))
  x_init <- runif(n_init, lob_x, upb_x)

  # initial update #####################################################
  for (i in 1:n_th){
    for (j in 1:n_init){
      # use negative since minimizing cycle time
      Y <- -QN_sim(parA, c(x_init[j], (2-x_init[j]), parS3), pr_m[i, , ], runlen, warmup)
      y <- mean(Y)

      # update parlist
      parlist_s <- update_s(parlist_s, i, x_init[j], y)
      parlist_d <- update_d(parlist_d, i, x_init[j], Y)

    }
  }

  # random choose parlist
  parlist_r <- parlist_s


  ###################################################################################
  # 4 Main procedure
  ###################################################################################
  results <- matrix(0, 0, 13)

  for (n_sim in 1:N){
    # calibration step ###################################################
    # summary
    lt_s <- apply(ltime(parlist_s$post_mu, parlist_s$post_va, target), 1, sum)
    idx_s <- which.max(lt_s)
    # detailed
    lt_d <- apply(ltime(parlist_d$post_mu, parlist_d$post_va, target), 1, sum)
    idx_d <- which.max(lt_d)
    # random
    idx_r <- sample((1:n_th), 1)


    # decision step ######################################################
    # summary
    x_s <- dec_s(parlist_s, idx_s, lob_x, upb_x)
    # detailed
    x_d <- dec_d(parlist_d, idx_d, lob_x, upb_x)
    # random
    x_r <- runif(1, lob_x, upb_x)


    # run simulation & update ############################################
    # summary
    y <- mean(-QN_sim(parA, c(x_s, (2-x_s), parS3), pr_m[idx_s, , ], runlen, warmup))
    parlist_s <- update_s(parlist_s, idx_s, x_s, y)
    # detailed
    Y <- -QN_sim(parA, c(x_d, (2-x_d), parS3), pr_m[idx_d, , ], runlen, warmup)
    parlist_d <- update_d(parlist_d, idx_d, x_d, Y)
    # random
    y_r <- mean(-QN_sim(parA, c(x_r, (2-x_r), parS3), pr_m[idx_r, , ], runlen, warmup))
    parlist_r <- update_s(parlist_r, idx_r, x_r, y_r)


    # final decision after current iteration #############################
    # summary
    res_s <- fdec(parlist_s, target, lob_x, upb_x)
    xhat_s <- res_s$x
    ihat_s <- res_s$th
    # detailed
    res_d <- fdec(parlist_d, target, lob_x, upb_x)
    xhat_d <- res_d$x
    ihat_d <- res_d$th
    # random
    res_r <- fdec(parlist_r, target, lob_x, upb_x)
    xhat_r <- res_r$x
    ihat_r <- res_r$th


    results <- rbind(results, c(n_sim, idx_s, x_s, ihat_s, xhat_s, 
                                idx_d, x_d, ihat_d, xhat_d,
                                idx_r, x_r, ihat_r, xhat_r))
  }

  # print(it)
  write.table(results, "results1.txt", sep=",", append=T, eol="\n",
            row.names=F, col.names=F);

}



