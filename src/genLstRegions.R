## Generate a list of ROIs --------------------
##  First use seed pixels to generate potential areas with signals, using genNewLstPixel(), 
## getting a binary map of potential pixels.
##  Then segment that binary map into consecutive ROIs, using regionGrow().
##  Output: lstRegions - A list of regions. Each element is a list of pixel coordinates within one region.
genLstRegions <- function(workZMap,thrZval){

    selMap<-workZMap*0
    availMap<-matrix(1,nrow=nrow(workZMap),ncol=ncol(workZMap)) # 1: available
    numSubRegion<-0
    mxZAvail<-max(workZMap[availMap==1])
    Zscore <- c()
    
    while(mxZAvail>thrZval*3/8){   # until the max z in available pixels cannot reject null hypothesis at 0.01 level
        numSubRegion<-numSubRegion+1
        cat(numSubRegion,'\n')
        seedPixel<-which(workZMap == mxZAvail, arr.ind = TRUE)  #seed pixel(s) has the max zScore among available pixels
        lstPixel<-list()
        for(i in 1:nrow(seedPixel)){   # now lstPixel is a list containing current batch of seeds pixels (their coordinates)
            lstPixel[[i]]<-seedPixel[i,]
        }
        res<-genNewLstPixel(lstPixel,workZMap,availMap)  # Expand the seeds to regions within available pixels
        newLstPixel<-res$newLstPixel
        Zscore<-c(Zscore,res$bestZ)
        for(i in 1:length(newLstPixel)){
            x<-newLstPixel[[i]][1]
            y<-newLstPixel[[i]][2]
            availMap[x,y]<-0
            if(res$bestZ>thrZval){
                selMap[x,y]<-1    # 1: already been used in some accepted region
            }
        } 
        mxZAvail<-max(workZMap[availMap==1])
        cat(mxZAvail,'\n')
    }
    
    if(0) {  # plot selMap
        idx <- which(workZMap>0, arr.ind = T)
        x0 <- min(idx[,1])-2
        x1 <- max(idx[,1])+2
        y0 <- min(idx[,2])-2
        y1 <- max(idx[,2])+2
        
        cat('bndr: ',paste(x0,x1,y0,y1,sep=' '),'\n')
        
        ## z-score map
        mat_z <- selMap
        mat_z <- mat_z[x0:x1,y0:y1]
        write.csv(mat_z,file=paste0(pth,'mat_region.csv'))
    }
    
    lstRegions<-list()
    availMap<-selMap  # availMap (as 1) means the selected non-noise areas (with some signals)
    m <- 0
    while(sum(availMap)>0){  # Each selected pixel is assigned to some ROI
        m <- m+1
        seedPixel<-which(workZMap == max(workZMap[availMap==1]), arr.ind = TRUE)
        lstPixel<-regionGrow(seedPixel,availMap)
        lstRegions[[m]]<-lstPixel
        for(i in 1:length(lstPixel)){
            x<-lstPixel[[i]][1]
            y<-lstPixel[[i]][2]
            availMap[x,y]<-0
        }
    }
    cat(Zscore,'\n')
    return(lstRegions)
}



