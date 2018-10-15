# Purpose: estimate VaR_beta of first M customers' times in M/M/1 queue
# Date : 1/30/2016           Author : Wei Xie

# Input:
#    arrival time ~ exp(parA), parA are arrival rates at k design points, (kx1)
#    service time ~ exp(parS), parS are service rates at k design points, (kx1)
#    nT - number of replications at design points
#    beta - consider beta-th percentile

# Ouput:
#    Y - percentiles of staying times (kxnT)

MM1_sim = function(parA,parS,runlen,warmup)
{
  runlength = runlen + warmup
  
  SE_times = matrix(0,runlength,2)  # arrival and depature times of first M customers
  temp = 0;                 # average number of customers in system
  arrN = 1;                 # record arrivals
  depN = 0;                 # record departures
  stat = 1;                 # number of customers in current system
      
  timer = rexp(1,rate=parA)  # first arrival 
  SE_times[arrN,1] = timer     # record the first customer's arrival time
  arrivalT = timer + rexp(1,rate=parA);
  serviceT = timer + rexp(1,rate=parS);

  while (depN < runlength){
        
    if (arrivalT < serviceT){   
      temp = temp + stat*(arrivalT-timer);
      timer = arrivalT;
      arrN = arrN+1;
      SE_times[arrN,1] = timer
      if (stat == 0){
        serviceT = timer+rexp(1,rate=parS);
      }
      stat = stat+1;
      if(arrN < runlength){
        arrivalT = timer+rexp(1,rate=parA);
      }else{
        arrivalT = Inf
      }
          
    }else{
      temp = temp + stat*(serviceT-timer);
      stat = stat-1;
      timer = serviceT;
      depN = depN+1;
      SE_times[depN,2] = timer
      if(stat>0){
        serviceT = timer+rexp(1,rate=parS);
      }else{
        serviceT = Inf;
      }
    }
  }                           # end of loop over runlength
#      temp2 = SE_times[,2]-SE_times[,1]    # staying times of first M customers
#      Y[m,j] = mean(temp2)
  temp2 = pmax(c(0,SE_times[(1:runlength-1),2])-SE_times[,1],0)
  Y = temp2[-(1:warmup)]
  
  return(Y);
}


