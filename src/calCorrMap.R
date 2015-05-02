## correlation between pixels and charX ------------------
calCorrMap_charX <- function(XData,lstPixels,shift,charX){
    corrMap<-matrix(0,nrow=dim(XData)[1],ncol=dim(XData)[2])
    numPixels<-length(lstPixels)
    vecLength <- dim(XData)[3]
    for(ii in 1:numPixels){
        x <- lstPixels[[ii]][1]    
        y <- lstPixels[[ii]][2]
        yCurve <- XData[x,y,]
        # shift0 <- shift[ii]
        shift0 <- 0
        if(shift0>0) {
            X <- charX[1:(vecLength-shift0)]
            Y <- yCurve[(shift0+1):vecLength]
        } else {
            X<-charX[(abs(shift0)+1):vecLength]
            Y<-yCurve[1:(vecLength-abs(shift0))]
        }
        corrMap[x,y] <- fast_cor(X,Y)
    }
    return(corrMap)
}

## correlation in neighbors of a matrix -------------
## Corr (x,Average wave) in each of the four directions
calCorrMap <- function(XData,lstPixels){
  corrMap<-matrix(0,nrow=dim(XData)[1],ncol=dim(XData)[2])
  numPixels<-length(lstPixels)
  for(i in 1:numPixels){
    x<-lstPixels[[i]][1]    
    y<-lstPixels[[i]][2]
    V0<-XData[x,y,]
    tempCorr<-rep(0,4)
    
    V1<-XData[x-1,y,]
    V2<-XData[x+1,y,]
    tempCorr[1]<-fast_cor(V0,(V1+V2)/2) # cor of avg, or avg of cor?
    
    V3<-XData[x,y-1,]
    V4<-XData[x,y+1,]
    tempCorr[2]<-fast_cor(V0,(V3+V4)/2)
    
    V5<-XData[x-1,y-1,]
    V6<-XData[x+1,y+1,]
    tempCorr[3]<-fast_cor(V0,(V5+V6)/2)
    
    V7<-XData[x+1,y-1,]
    V8<-XData[x-1,y+1,]
    tempCorr[4]<-fast_cor(V0,(V7+V8)/2)
    
    corrMap[x,y]<-max(tempCorr,na.rm = TRUE)
    # if(corrMap[x,y] < -1) corrMap[x,y]<-0   #if tempoCorr is empty
  }
  
  # corrMap[is.na(corrMap)]<-0  # let the edge points of corrMap to be 0
  return(corrMap)
}

## correlation in neighbors of a matrix -------------
## use absolute position vector as input
calCorrMap1 <- function(XData,vecPixels){
    # Corr (x,Average wave) in each of the four directions
    corrMap<-matrix(0,nrow=dim(XData)[1],ncol=dim(XData)[2])
    numPixels<-length(vecPixels)
    nx <- dim(XData)[1]
    nxp1 <- nx + 1
    nxm1 <- nx - 1
    ny <- dim(XData)[2]
    if(nx!=ny) {
        stop('oops\n')
    }
    NN <- dim(XData)[1] * dim(XData)[2]
    rg <- 0:(dim(XData)[3]-1)
    for(ii in 1:numPixels){
        xx <- vecPixels[ii]
        idx <- xx + NN*rg
        V0<-XData[idx]
        tempCorr<-rep(0,4)
        
        V1<-XData[idx-1]
        V2<-XData[idx+1]
        tempCorr[1]<-fast_cor(V0,(V1+V2)/2) # cor of avg, or avg of cor?
        
        V3<-XData[idx-nx]
        V4<-XData[idx+nx]
        tempCorr[2]<-fast_cor(V0,(V3+V4)/2)
        
        V5<-XData[idx-nxp1]
        V6<-XData[idx+nxp1]
        tempCorr[3]<-fast_cor(V0,(V5+V6)/2)
        
        V7<-XData[idx-nxm1]
        V8<-XData[idx+nxm1]
        tempCorr[4]<-fast_cor(V0,(V7+V8)/2)
        
        corrMap[xx]<-mean(tempCorr,na.rm = TRUE)
    }
    return(corrMap)
}

