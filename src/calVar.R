calVar <- function(totN, ordSelN){
  ordSelN<-sort(ordSelN)
  tempVar<-0
  Finv<--qnorm(c(1:(totN+1))/(totN+2))
  difFJ<-Finv[1:totN]-Finv[(1:totN)+1]
  
  #show(difFJ)
  #show(sum(difFJ))
  for(i in 1:length(ordSelN)){
    difFI<-difFJ[ordSelN[i]]
    tempVar<-tempVar+2*sum(ordSelN[i]/(totN+1)*(1-ordSelN[i:length(ordSelN)]/(totN+1))
                           *difFI*difFJ[ordSelN[i:length(ordSelN)]])
  }
  return(tempVar*totN)
}
