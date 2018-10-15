#require(actuar)
require(SpatialExtremes)

# Service: Generalized Pareto distribution with loc, scale, shape

MP1_sim = function(parA,parS_lc,parS_sc,parS_sh,runlen,warmup) # Service: Generalized Pareto distribution
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
  serviceT = timer + rgpd(1, loc=parS_lc, scale=parS_sc, shape=parS_sh); # rgpd(n, loc, scale, shape)

  while (depN < runlength){
        
    if (arrivalT < serviceT){   
      temp = temp + stat*(arrivalT-timer);
      timer = arrivalT;
      arrN = arrN+1;
      SE_times[arrN,1] = timer
      if (stat == 0){
        serviceT = timer+rgpd(1, loc=parS_lc, scale=parS_sc, shape=parS_sh);
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
        serviceT = timer+rgpd(1, loc=parS_lc, scale=parS_sc, shape=parS_sh);
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


