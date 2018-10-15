# Purpose: output first M customers' times in queueing network
# Date : 8/18/2017           Author : Bo Wang

# Input:
#    arrival time ~ exp(parAi), parAi are arrival rates at i queue
#    service time ~ exp(parSi), parSi are service rates at i queue


# Ouput:
#    Y - staying times 

QN_sim = function(parA,parS,pr,runlen,warmup)
{
  #parA arrival rate at 1
  #parS = c(parS1,parS2,parS3) service rate at 1,2,3
  #pr = probability matrix

  runlength = runlen + warmup

  SE_times <- matrix(0,0,2) # record arrival and departure times of customers
  s1_q <- s2_q <- s3_q <- c() # queue for stations 1,2,3

  arrN <- 1 # arrival number (also an index)
  timer <- rexp(1,rate=parA) # first arrival
  SE_times <- rbind(SE_times,c(timer,0)) # record first arrival time
  s1_q <- c(s1_q,arrN) # append first customer into queue 1

  T1_arr <- timer + rexp(1,rate=parA); # next arrival at station 1
  T1_out <- timer + rexp(1,rate=parS[1]); # next out at station 1
  T2_out <- Inf; # next out at station 2
  T3_out <- Inf; # next out at station 3
  next_one <- c( T1_arr, T1_out, T2_out, T3_out);

  while (arrN < 1000000){

    if (which.min(next_one) == 1){ #(C1) arrival at 1
      arrN <- arrN + 1
      timer <- next_one[1]
      SE_times <- rbind(SE_times,c(timer,0))

      if (length(s1_q) == 0){
        T1_out <- timer+rexp(1,rate=parS[1])
        next_one[2] <- T1_out
      }
      s1_q <- c(s1_q,arrN)

      T1_arr <- timer + rexp(1,rate=parA)
      next_one[1] <- T1_arr
    
    }else if (which.min(next_one) == 2){ #(C2) out of 1
      timer <- next_one[2]

      # determine queue 2 or queue 3 to go
      if (runif(1) < pr[1,2]){
        if (length(s2_q) == 0){
          T2_out <- timer+rexp(1,rate=parS[2])
          next_one[3] <- T2_out
        }
        s2_q <- c(s2_q,s1_q[1])
      }else{
        if (length(s3_q) == 0){
          T3_out <- timer+rexp(1,rate=parS[3])
          next_one[4] <- T3_out
        }
        s3_q <- c(s3_q,s1_q[1])
      }
      s1_q <- s1_q[-1]

      if (length(s1_q) > 0){
        T1_out <- timer+rexp(1,rate=parS[1])
      }else{
        T1_out <- Inf
      }
      next_one[2] <- T1_out

    }else if (which.min(next_one) == 3){ #(C3) out of 2
      timer <- next_one[3]

      # determine go back to 1 or depart
      if (runif(1) < pr[2,1]){
        if (length(s1_q) == 0){
          T1_out <- timer + rexp(1,rate=parS[1])
          next_one[2] <- T1_out
        }
        s1_q <- c(s1_q,s2_q[1])
      }else{
        SE_times[s2_q[1],2] <- timer
      }
      s2_q <- s2_q[-1]

      if (length(s2_q) > 0){
        T2_out <- timer + rexp(1,rate=parS[2])
      }else{
        T2_out <- Inf
      }
      next_one[3] <- T2_out
    
    }else if (which.min(next_one) == 4){ #(C4) out of 3
      timer <- next_one[4]

      # determine go back to 1 or depart
      if (runif(1) < pr[3,1]){
        if (length(s1_q) == 0){
          T1_out <- timer + rexp(1,rate=parS[1])
          next_one[2] <- T1_out
        }
        s1_q <- c(s1_q,s3_q[1])
      }else{
        SE_times[s3_q[1],2] <- timer
      }
      s3_q <- s3_q[-1]

      if (length(s3_q) > 0){
        T3_out <- timer + rexp(1,rate=parS[3])
      }else{
        T3_out <- Inf
      }
      next_one[4] <- T3_out

    }

    # break loop if first (runlength) customers have departed
    if (length(SE_times[,2]) >= runlength && 
      min(SE_times[(1:runlength),2]) > 0){
      break;
    }
  
  }
  # compute cycle times 
  cycleT <- SE_times[(1:runlength),2] - SE_times[(1:runlength),1]
  Y <- cycleT[-(1:warmup)]
  
  return(Y)
}




