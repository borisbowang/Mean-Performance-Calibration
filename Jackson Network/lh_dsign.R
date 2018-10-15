require(lhs)

set.seed(10) #set seed

# service rate
parS <- maximinLHS(30,3)*1.4+0.8
write.table(parS, "parS_d1.txt")




