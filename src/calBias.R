calBias <- function(totN, ordSelN){
  bias<-0
  for(i in 1:length(ordSelN)){
    bias<-bias-qnorm(ordSelN[i]/(totN+1))
  }
  return(bias)
}
