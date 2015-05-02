checkEnclose1 <- function(pixel,occMap){
  # Check whether current pixel is "enclosed"/surrounded by pixels in list (what occMap indicates)
  # return enclosed: 1-Y, 0-N
  x<-pixel[1]
  y<-pixel[2]
  enclosed <- 0
  
  y1 <- 512*(y-1)
  y1l <- y1 - 512
  y1r <- y1 + 512
  
  gp1 <- c(y1+x+1,y1+x-1,y1l+x,y1r+x)
  numEnclosedPoints <- sum(occMap[gp1]) 
  if(numEnclosedPoints>3) enclosed<-1                              #???  Why 3?
  gp2 <- c(y1l+x+1,y1l+x-1,y1r+x+1,y1r+x-1)
  numEnclosedPoints <- numEnclosedPoints + sum(occMap[gp2]) 
  if(numEnclosedPoints>4) enclosed<-1                              #??? Why totally 4?

  return(enclosed)
}

checkEnclose <- function(pixel,occMap){
    # Check whether current pixel is "enclosed"/surrounded by pixels in list (what occMap indicates)
    # return enclosed: 1-Y, 0-N
    x<-pixel[1]
    y<-pixel[2]
    
    enclosed<-0
    numEnclosedPoints<-0
    if(occMap[x-1,y]==1) numEnclosedPoints<-numEnclosedPoints+1
    if(occMap[x+1,y]==1) numEnclosedPoints<-numEnclosedPoints+1
    if(occMap[x,y-1]==1) numEnclosedPoints<-numEnclosedPoints+1
    if(occMap[x,y+1]==1) numEnclosedPoints<-numEnclosedPoints+1
    if(numEnclosedPoints>3) enclosed<-1                              #???  Why 3?
    
    if(occMap[x-1,y-1]==1) numEnclosedPoints<-numEnclosedPoints+1
    if(occMap[x-1,y+1]==1) numEnclosedPoints<-numEnclosedPoints+1
    if(occMap[x+1,y-1]==1) numEnclosedPoints<-numEnclosedPoints+1
    if(occMap[x+1,y+1]==1) numEnclosedPoints<-numEnclosedPoints+1
    
    if(numEnclosedPoints>4) enclosed<-1                              #??? Why totally 4?
    
    return(enclosed)
}


