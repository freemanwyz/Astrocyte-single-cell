

genNewLstPixel <- function(lstPixel,workZMap,availMap){
  ### Enlarge the region from one/several seed pixel(s)
  #   In each iteration, get all adjacent neighbor pixels of current region, calculate accumulated z's,
  #   find out the set of neighbors which lead to the highest mu*sqrt(n), and add them into the list of pixels
  #   Input:  lstPixel - A list of seed pixels
  #           workZMap - A map of z scores
  #           availMap - A map indicating available pixels
  #   Output: res
  #             res$newLstPixel - The union of seed pixels and found neighbors
  #             res$bestZ - The highest z score (???) of this region
  newLstPixel <- lstPixel
  #m<-0
  numNewPoints<-0
  numOldNeigbors<-0
  bestZ<-workZMap[lstPixel[[1]][1],lstPixel[[1]][2]]
  while(TRUE){
    #m<-m+1                                          # This m is surely unused before it's changed by another command. Delete it???
    lstNeighbor<-findNeighbor(newLstPixel,availMap)  # Find all direct neighbors(who're in availMap) of current pixels in newLstPixel 
    #show(length(newLstPixel)+length(lstNeighbor))
    if(length(lstNeighbor)==0) break  # if no more neighbors can be found, done
    
    ZSaN<-rep(0,length(newLstPixel)+length(lstNeighbor))
    for(i in 1:length(newLstPixel)){   # the first part of ZSaN is consist of z of pixels in newLstPixels
      x<-newLstPixel[[i]][1]
      y<-newLstPixel[[i]][2]
      ZSaN[i]<-workZMap[x,y]
    }
    
    ## Sort lstNeighbor by z's, get the ordered z's as stZs (descending order)
    Zs<-c(1:length(lstNeighbor))
    for(i in 1:length(lstNeighbor)){
      x<-lstNeighbor[[i]][1]
      y<-lstNeighbor[[i]][2]
      Zs[i]<-workZMap[x,y]
    }
    rankZ<-order(-Zs)
    lstNeighbor<-lstNeighbor[rankZ]  #list of current neighbors being ranked by z-score
    stZs<-Zs[rankZ]  # corresponding ranked z-scores of these neighbors
 
    ############################################################################
    ## Set numEquivPix                                                              ???
    numEquivPix<-rep(0,1+length(lstNeighbor))
    numSigPix<-length(newLstPixel)
    for(i in 1:length(numEquivPix)){
      if((numSigPix+i-1)==1){
        numEquivPix[i]<-1
      }
      else{
        numEquivPix[i]<-1.5*(numSigPix+i-1)
      }
      numEquivPix[i]<-numSigPix+i-1           #???                                              #???  
      #Perhaps this line should be deleted? Otherwise the above "if" is meaningless!
    }
    
    ## Let smallZs[1] be the sum of z-scores of all pixels in newLstPixel 
    smAllZs<-rep(0,1+length(lstNeighbor))
    for(i in 1:length(newLstPixel)){ 
      x<-newLstPixel[[i]][1]
      y<-newLstPixel[[i]][2] 
      smAllZs[1]<-smAllZs[1]+workZMap[x,y]
    }
    
    occMap<-workZMap*0
    for(i in 1:length(newLstPixel)){  # occMap: indicator matrix of pixels in newLstPixel
      x<-newLstPixel[[i]][1]
      y<-newLstPixel[[i]][2] 
      occMap[x,y]<-1
    }
    
    ## Rearrange the neighbors in lstNeighbor, putting the rearranged ones in updLstNeighbor         #But why???
    ## Put the corresponding z's to ZSaN, in the same order as in updLstNeighbor
    ## Add elements corresponding to this batch of neighbors to smAllZs                              #including 0s???
    ## Set occMap(pixels in this batch of neighbors) <- 1
    flagNeighbor<-rep(0,length(lstNeighbor))
    updLstNeighbor<-lstNeighbor   # What about just use 0s so that we can see the program more clearly      ???
    m<-0                             # Wait, wait... We set m to be 0 every loop?                               ???
    for(i in 1:length(lstNeighbor)){  # for each currently potential neighbor (in lstNeighbor) i
      # smallZs is an accumulating series: smallZs[k]=smallZs[k-1]+z[k] (k is an index in updLstNeighbor)
      if(flagNeighbor[i]==0){
        m<-m+1
        updLstNeighbor[m]<-lstNeighbor[i]  # if flag is 0, put this potential neighbor to updLstNeighbor
        tempZ<-smAllZs[m]+stZs[i]  # sum of z's of pixels in newLstPixel and updLstNeighbor
        flagNeighbor[i]<-1
        x<-lstNeighbor[[i]][1]
        y<-lstNeighbor[[i]][2] 
        occMap[x,y]<-1
        ZSaN[numSigPix+m]<-workZMap[x,y]    # the second part of ZSaN are z of current neighbors (updLstNeighbor)
        
## What is the purpose of the following for loop ??? --> find the enclosed pixels after we take into account the above pixel
        for(j in (i+1):length(lstNeighbor)){
          if(i+1>length(lstNeighbor)) next  # How about "break"  ???  The condition is the same for any j
          
          if(flagNeighbor[j]==0){  # if flag of j is 0 (j=i+1:length(potential neighbor))
            pixel<-lstNeighbor[[j]]
            # browser()
            enclosed<-checkEnclose(pixel,occMap)  # check whether j is enclosed by occed pixels  ??? (criteria???)
            if(enclosed==1){  # if j IS enclosed
              m<-m+1
              updLstNeighbor[m]<-lstNeighbor[j]
              tempZ<-tempZ+stZs[j]
              flagNeighbor[j]<-1
              x<-lstNeighbor[[j]][1]
              y<-lstNeighbor[[j]][2] 
              occMap[x,y]<-1
              ZSaN[numSigPix+m]<-workZMap[x,y]
            }
          }
        }
        smAllZs[m+1]<-tempZ 
        # So... if there IS some enclosed neighbors, there will be some 0s in smAllZs, correct? --> Yes
      }
    }

    numNewNeigbors<-sum(smAllZs[2:length(smAllZs)]>0) #Only interested in those non-zero elements in smAllZs? Don't we include enclosed ones???
    ordZSaN<-order(-ZSaN)
    biasAllZs<-rep(0,length(smAllZs))    # What do you mean by "variance" and "bias" here     ???
    varAllZs<-rep(0,length(smAllZs))
    totN<-length(newLstPixel)+length(lstNeighbor)
    facVar<-totN/calVar(totN,c(1:totN))              # ???
    for(i in 1:length(biasAllZs)){                   # ??? This for loop has no function now, right?
      selN<-length(newLstPixel)+(i-1)
      #biasAllZs[i]<-calBias(totN,ordZSaN[1:selN])
      #varAllZs[i]<-calVar(totN,ordZSaN[1:selN])*facVar
    }
    
    #neighborZ<-(smAllZs-biasAllZs)/sqrt(varAllZs)
    neighborZ<-smAllZs/sqrt(numEquivPix)            # Why not record the number of pixels when including in new neighbors???
                                          # 'Cause the number of pixels can be varying based on the "enclosed" condition, can't it???
    
    numOldNeigbors<-numNewNeigbors
    
    ############################################################################
    #show(workZMap[newLstPixel[[1]][1],newLstPixel[[1]][2]])
    #show(stZs)
    #show(smAllZs)
    #show(biasAllZs)
    #show(varAllZs)
    #show(neighborZ)
    #readline("pause; enter to continue")
    numNewPoints<-order(-neighborZ+0)[1]   # find the index of the max mu*sqrt(n) in neighborZ/smAllZs
    
    selN<-length(newLstPixel)+(numNewPoints-1)
    biasBestZ<-calBias(totN,ordZSaN[1:selN])    #   ???
    varBestZ<-calVar(totN,ordZSaN[1:selN])*facVar    #   ???
    tempBestZ<-(smAllZs[numNewPoints]-biasBestZ)/sqrt(varBestZ)   #   ???
    #show(tempBestZ)
    #show(neighborZ[numNewPoints])

    if(numNewPoints==1){     # if the max mu*sqrt(n) appears at smAllZs[1], say, not including any neighbor
      break
    }
    else{
      if(tempBestZ<=bestZ){
        tempnumNewPoints<-numNewPoints
        for(i in numNewPoints:2){
          selN<-length(newLstPixel)+(i-1)
          revalBestZ<-smAllZs[i]-calBias(totN,ordZSaN[1:selN])                 #  ???
          revalBestZ<-revalBestZ/sqrt(calVar(totN,ordZSaN[1:selN])*facVar)     #  ???
          revalBestZ<revalBestZ+qnorm(1/(numNewPoints-i+2))                    #  ???
          #show("revalBestZ")
          #show(revalBestZ)
          if(bestZ < revalBestZ){
            tempBestZ<-revalBestZ
            tempnumNewPoints<-i
            numNewPoints<-tempnumNewPoints
            break
          }
        }
      }
      
      if(tempBestZ<=bestZ){
        break
      }
      bestZ<-tempBestZ
      currLength<-length(newLstPixel)
      for(i in 2:numNewPoints){
        newLstPixel[currLength+i-1]<-updLstNeighbor[i-1]
      }
    }  
  }
  
  res<-NULL
  res$newLstPixel<-newLstPixel
  res$bestZ<-bestZ
  return(res)
}



