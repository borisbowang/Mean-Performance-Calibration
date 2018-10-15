rm(list=ls())
source("QN_sim.R")
source("Update.R")

# require(ggplot2)
# require(reshape2)

set.seed(10) #set seed
#########################################################################
#1 setting M alternatives (adjust here)
#########################################################################
# runlength
runlen <- 50
# warmup
warmup <- 1000

# simulation budget
N <- 500
# micro replications     
IT <- 150
   

# arrival rate = 0.5
parA <- 0.5
# probability
pr <- rbind(c(0, 0.7, 0.3, 0),
      c(0.3, 0, 0, 0.7),
      c(0.2, 0, 0, 0.8))

arr1 <- parA/(1-pr[1,2]*pr[2,1]-pr[1,3]*pr[3,1])
arr2 <- arr1*pr[1,2]
arr3 <- arr1*pr[1,3]


# service rate (generate by Latin-hyper-cube min-max design)
parS <- as.matrix(read.table("parS_d1.txt"))

# num of alternatives   
M <- nrow(parS)


# true mean cycle time and utilization at 1 (M/M/1 case)
m_cyc <- rep(0,M)
# utl <- arr1/parS1

for (j in 1:M){
  m_cyc[j] <- sum( c(arr1,arr2,arr3)/(parS[j,]-c(arr1,arr2,arr3)) )/parA
}

# parS1
# utl
# sort(m_cyc)


# target mean cycle time
target <- 2


true.mse <- (m_cyc-target)^2

trueid <- which.min(true.mse)


# num of initial runs for each alternative
init_r <- 5

# repeat to estimate empirical PCS and mean OC
results <- matrix(0, 0, 7)

# # record the detailed results
# re_s <- matrix(0, 0, 2*M+3)
# re_kn <- matrix(0, 0, 2*M+3)
# re_ar <- matrix(0, 0, 2*M+3)
# re_sw <- matrix(0, 0, 2*M+3)

#########################################################################
#2 setting prior
#########################################################################

q <- tau <- rep(1,M)
al <- a <- rep(3/2,M)

# true (long-run estimation) Rr and AR(1) estimated Rr
Rr <- array(0,dim=c(runlen,runlen,M))
for (j in 1:M){
  yt <- as.ts( QN_sim(parA,parS[j,],pr,5000,warmup) )
  sacf <- as.vector(acf(yt,lag.max=(runlen-1),type="correlation",plot=F)$acf)
  R_j <- toeplitz(sacf)
  Rr[,,j] <- solve(R_j)
}



for (it in 1:IT)
{
  # trupar <- rbind(muv, r0v, phiv)
  theta_d <- theta_s <- rnorm(M)
  # q <- tau <- rep(1,M)
  # al <- a <- rep(3/2,M)
  be <- b <- runif(M, 200, 500)
  var_d <- var_s <- b/(tau*(a-1))

  # # true Rr
  # Rr <- array(0,dim=c(runlen,runlen,M))
  # for (j in 1:M){ Rr[,,j] <- RR(phiv[j],runlen) }

  # summary approach
  par_s0 <- list(theta_s=theta_s,tau=tau,a=a,b=b,var_s=var_s)

  # detailed approach with known Rr
  par_d_kn0 <- list(theta_d=theta_d,q=q,al=al,be=be,var_d=var_d)

  # detailed approach with known AR(1) Rr
  par_d_ar0 <- list(theta_d=theta_d,q=q,al=al,be=be,var_d=var_d)

  # detailed approach with stepwise rho estimation
  Ydat <- list()
  for (j in 1:M){ Ydat[[j]] <- matrix(0,0,runlen) }
  par_d_sw0 <- list(theta_d=theta_d,
                    q=q,
                    Ydat=Ydat,
                    al=al,
                    be=be,
                    var_d=var_d)


#########################################################################
#3 Numerical Study
#########################################################################


  par_s <- par_s0
  par_d_kn <- par_d_kn0
  par_d_sw <- par_d_sw0
  ##################################################
  # initial runs for each alternative
  ##################################################
  Ydat <- par_d_sw$Ydat
  if (init_r > 0){
    for (j in 1:M){
      for (i in 1:init_r){
        Y <- QN_sim(parA,parS[j,],pr,runlen,warmup)
        ybar <- mean(Y)

        par_s <- updates(par_s, j, ybar)
        par_d_kn <- updated_kn(par_d_kn, j, Rr[,,j], Y)

        Ydat[[j]] <- rbind(Ydat[[j]],Y)
      }

      Rr_j <- RR( rho_hat(Ydat[[j]]), runlen)

      for (i in 1:init_r){
        par_d_sw <- updated_kn(par_d_sw, j, Rr_j, Ydat[[j]][i,])
      }
    }
    par_d_sw$Ydat <- Ydat
  }

  ##################################################
  # adaptive learning
  ##################################################
  for (i in 1:N){
    #choose the alternative to run simulation--use local time
    #summary
    lt_s <- ltime(par_s$theta_s, par_s$var_s, target)
    idx_s <- r_select(which(lt_s == max(lt_s)))
    #detailed with known Rr
    lt_kn <- ltime(par_d_kn$theta_d, par_d_kn$var_d, target)
    idx_d_kn <- r_select(which(lt_kn == max(lt_kn)))
    #detailed with stepwise rho
    lt_sw <- ltime(par_d_sw$theta_d, par_d_sw$var_d, target)
    idx_d_sw <- r_select(which(lt_kn == max(lt_kn)))


    #generate data & update parameters
    #detailed with known Rr
    Y_kn <- QN_sim(parA,parS[idx_d_kn,],pr,runlen,warmup)
    par_d_kn <- updated_kn(par_d_kn, idx_d_kn, Rr[,,idx_d_kn], Y_kn)

    #detailed with stepwise rho
    if (idx_d_sw == idx_d_kn){
      Y_sw <- Y_kn
    }else{
      Y_sw <- QN_sim(parA,parS[idx_d_sw,],pr,runlen,warmup)
    }
    par_d_sw <- updated_sw(par_d_sw, idx_d_sw, Y_sw)

    #summary stat
    if (idx_s == idx_d_kn){
      ybar <- mean(Y_kn)
    }else if (idx_s == idx_d_sw){
      ybar <- mean(Y_sw)
    }else{
      ybar <- mean(QN_sim(parA,parS[idx_s,],pr,runlen,warmup))
    }
    par_s <- updates(par_s, idx_s, ybar)


    #choose the best alternative & record
    mse.s <- (par_s$theta_s-target)^2 + par_s$var_s
    mse.d_kn <- (par_d_kn$theta_d-target)^2 + par_d_kn$var_d
    mse.d_sw <- (par_d_sw$theta_d-target)^2 + par_d_sw$var_d

    results <- rbind(results, 
      c(i, r_select(which(mse.s==min(mse.s))), r_select(which(mse.d_kn==min(mse.d_kn))),
       r_select(which(mse.d_sw==min(mse.d_sw))), idx_s, idx_d_kn, idx_d_sw))
    
    # # record the detailed results
    # re_s <- rbind(re_s,
    #   c(i, idx_s, rbind(par_s$theta_s, par_s$var_s), r_select(which(mse.s==min(mse.s)))))
    # re_kn <- rbind(re_kn,
    #   c(i, idx_d_kn, rbind(par_d_kn$theta_d, par_d_kn$var_d), r_select(which(mse.d_kn==min(mse.d_kn)))))
    # re_ar <- rbind(re_ar,
    #   c(i, idx_d_ar, rbind(par_d_ar$theta_d, par_d_ar$var_d), r_select(which(mse.d_ar==min(mse.d_ar)))))
    # re_sw <- rbind(re_sw,
    #   c(i, idx_d_sw, rbind(par_d_sw$theta_d, par_d_sw$var_d), r_select(which(mse.d_sw==min(mse.d_sw)))))

  }

  print (it)
}


################################################################
#2 plot results
################################################################

results <- as.data.frame(results)
names(results) <- c("Nstep", "best.s", "best.d_kn", "best.d_sw",
 "design.s", "design.d_kn", "design.d_sw")

results$cs.summary <- 1*(results$best.s==trueid)
results$cs.detail_kn <- 1*(results$best.d_kn==trueid)
results$cs.detail_sw <- 1*(results$best.d_sw==trueid)

results$cost.summary <- true.mse[results$best.s]-true.mse[trueid]
results$cost.detail_kn <- true.mse[results$best.d_kn]-true.mse[trueid]
results$cost.detail_sw <- true.mse[results$best.d_sw]-true.mse[trueid]


# pcs <- aggregate(cbind(cs.summary, cs.detail_kn, cs.detail_sw)~Nstep, 
#                             data=results, mean)
# pcs <- melt(pcs, id="Nstep")

# p <- ggplot(pcs, aes(Nstep, value, group=variable,   
#   color=variable))+
#   geom_line() + theme_bw()

# ggsave(paste("M_",M,"pcs.pdf", sep=""))  


# pcost <- aggregate(cbind(cost.summary, cost.detail_kn, cost.detail_sw)~Nstep, 
#                             data=results, mean)
# pcost <- melt(pcost, id="Nstep")
# p1 <- ggplot(pcost, aes(Nstep, value, group=variable,   
#   color=variable))+
#   geom_line() + theme_bw()

# ggsave(paste("M_",M,"pcost.pdf", sep=""))  


# # the detailed results
# re_s <- as.data.frame(re_s)
# re_kn <- as.data.frame(re_kn)
# re_ar <- as.data.frame(re_ar)
# re_sw <- as.data.frame(re_sw)
# names(re_s) <- names(re_kn) <- names(re_ar) <- names(re_sw) <- c("Nstep", "design", 
#   c(rbind(paste("mu_",1:M), paste("var_",1:M))),
#  "best")



# write.table(pcost,"pcost.txt")
# write.table(pcs,"pcs.txt")
#write.table(results,"results2.txt",append=T,row.names=F)

write.table(results,"results2.txt",sep=",",append=T,eol="\n",
            row.names=F,col.names=F)


# #summary 
# pcost[c(50,100,200,500),3]
# #detail kn
# pcost[c(550,600,700,1000),3]
# #detail sw
# pcost[c(1050,1100,1200,1500),3]



