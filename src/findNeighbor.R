findNeighbor <- function(lstPixel,availMap){
  # For any pixel in current lstPixel list, check its 8 neighbor pixels
  # If a neighbor is in available map, add it to neighbor list
  # Return the list of neighbors
  
  lstNeighbor<-list()
  availMap[1,]<-0  # the borders should not be included in the available map
  availMap[nrow(availMap),]<-0
  availMap[,1]<-0
  availMap[,ncol(availMap)]<-0
  
  for(i in 1:length(lstPixel)){  # pixels in lstPixel have been used already, shouldn't be included in availabel map
    x<-lstPixel[[i]][1]
    y<-lstPixel[[i]][2]
    availMap[x,y]<-0
  }
  
  for(i in 1:length(lstPixel)){
    # For any pixel in current lstPixel list, check its 8 neighbor pixels
    # If the neighbor is in available map, add it to neighbor list and remove it from available (only for this function) map
    x<-lstPixel[[i]][1]
    y<-lstPixel[[i]][2]
  
    tempX<-x-1
    tempY<-y-1
    if(availMap[tempX,tempY]==1){
      curLength<-length(lstNeighbor)
      lstNeighbor[[curLength+1]]<-c(tempX,tempY)
      availMap[tempX,tempY] <- 0
    }
    
    tempX<-x-1
    tempY<-y
    if(availMap[tempX,tempY]==1){
      curLength<-length(lstNeighbor)
      lstNeighbor[[curLength+1]]<-c(tempX,tempY)
      availMap[tempX,tempY] <- 0
    }
    
    tempX<-x-1
    tempY<-y+1
    if(availMap[tempX,tempY]==1){
      curLength<-length(lstNeighbor)
      lstNeighbor[[curLength+1]]<-c(tempX,tempY)
      availMap[tempX,tempY] <- 0
    }
    
    tempX<- x
    tempY<- y-1
    if(availMap[tempX,tempY]==1){
      curLength<-length(lstNeighbor)
      lstNeighbor[[curLength+1]]<-c(tempX,tempY)
      availMap[tempX,tempY] <- 0
    }
    
    tempX<-x
    tempY<-y+1
    if(availMap[tempX,tempY]==1){
      curLength<-length(lstNeighbor)
      lstNeighbor[[curLength+1]]<-c(tempX,tempY)
      availMap[tempX,tempY] <- 0
    }
    
    tempX<-x+1
    tempY<-y-1
    if(availMap[tempX,tempY]==1){
      curLength<-length(lstNeighbor)
      lstNeighbor[[curLength+1]]<-c(tempX,tempY)
      availMap[tempX,tempY] <- 0
    }
    
    tempX<-x+1
    tempY<-y
    if(availMap[tempX,tempY]==1){
      curLength<-length(lstNeighbor)
      lstNeighbor[[curLength+1]]<-c(tempX,tempY)
      availMap[tempX,tempY] <- 0
    }
    
    tempX<-x+1
    tempY<-y+1
    if(availMap[tempX,tempY]==1){
      curLength<-length(lstNeighbor)
      lstNeighbor[[curLength+1]]<-c(tempX,tempY)
      availMap[tempX,tempY] <- 0
    }    
  }
  return(lstNeighbor)
}
