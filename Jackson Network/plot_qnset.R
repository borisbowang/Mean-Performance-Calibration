rm(list=ls())
require(lhs)

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


# service rate
set.seed(10)
parS <- maximinLHS(30,3)*1.4 + 0.8

# num of alternatives   
M <- nrow(parS)


# true mean cycle time and utilization at 1 (M/M/1 case)
m_cyc <- rep(0,M)
# utl <- arr1/parS1

for (j in 1:M){
  m_cyc[j] <- sum( c(arr1,arr2,arr3)/(parS[j,]-c(arr1,arr2,arr3)) )/parA
}

epc <- sort(m_cyc)

pdf("scatterplots.pdf", width=6, height=4.5)
plot(epc,type = 'o',ylim = c(1,7),xlab = "Calibration Settings Index",
     ylab = "Expected Cycle Time", main = "Jackson Network Settings")
lines(c(-1,30),c(2,2),lty=2,col="blue")
lines(c(-1,30),c(3,3),lty=2,col="blue")
lines(c(-1,30),c(4,4),lty=2,col="red")
lines(c(-1,30),c(5,5),lty=2,col="blue")
lines(c(-1,30),c(6,6),lty=2,col="blue")
legend("topleft",c("Candidate Calibration Settings"),
       pch=c(1))
dev.off()






