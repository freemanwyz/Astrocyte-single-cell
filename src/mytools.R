## region grow ----------
my_region_grow <- function(dat,seed,thr=0.75) {
    neib0 <- c(1,0,-1,0,0,1,0,-1,1,1,1,-1,-1,1,-1,-1)
    neib0 <- matrix(neib0,nrow = 8,ncol = 2,byrow = T)
    Nx <- dim(dat)[2]
    lst <- matrix(seed,nrow=1,ncol=2)
    mat_avail <- dat*0+1
    mat_avail[lst] <- 0
    thry <- dat[lst]*thr
    ii <- 1
    while(ii<=dim(lst)[1]) {
        ## neighbor within the boundary
        neib <- t(t(neib0) + as.vector(lst[ii,]))
        ii <- ii + 1
        neib <- my_within_boundary(neib,Nx)
        if(dim(neib)[1]==0) next
        ## neighbor not searched yet
        idx <- mat_avail[neib]>0
        neib <- neib[idx,,drop=F]
        if(dim(neib)[1]==0) next
        ## neighbor higher than threshold
        mat_avail[neib] <- 0
        idx <- dat[neib] >= thry
        neib <- neib[idx,,drop=F]
        if(dim(neib)[1]==0) next
        ## new list of pixels
        lst <- rbind(lst,neib)
    }
    return(lst)
}

my_region_grow_n <- function(dat,lst,n) {
    for(ii in seq(n)) {
        lst <- my_region_grow_one(dat,lst)
    }
    return(lst)
}

my_region_grow_one <- function(dat,lst) {
    neib0 <- c(1,0,-1,0,0,1,0,-1,1,1,1,-1,-1,1,-1,-1)
    neib0 <- matrix(neib0,nrow = 8,ncol = 2,byrow = T)
    Nx <- dim(dat)[2]
    mat_avail <- dat*0+1
    mat_avail[lst] <- 0
    for(ii in seq(dim(lst)[1])) {
        ## neighbor within the boundary
        neib <- t(t(neib0) + as.vector(lst[ii,]))
        ii <- ii + 1
        neib <- my_within_boundary(neib,Nx)
        if(dim(neib)[1]==0) next
        ## neighbor not searched yet
        idx <- mat_avail[neib]>0
        neib <- neib[idx,,drop=F]
        if(dim(neib)[1]==0) next
        mat_avail[neib] <- 0
        ## new list of pixels
        lst <- rbind(lst,neib)
    }
    return(lst)
}

## load data ---------
my_load_dat <- function() {
    ## neuron
    img <- readImage(datfp$ap,type = 'tiff',all = T)
    img <- channel(img,'green')
    img <- img@.Data
    f_max <- apply(img,c(1,2),max)
    f_bar <- apply(img,c(1,2),mean)
    img0 <- f_max-f_bar
    img0s <- img0/max(img0)
    ## synapse
    img1 <- readImage(datfp$maxz_synapsin,type = 'tiff',all = T)
    img1 <- channel(img1,'gray')
    img1s <- img1/max(img1)
    img2 <- readImage(datfp$maxz_PSD95,type = 'tiff',all = T)
    img2 <- channel(img2,'gray')
    img2s <- img2/max(img2)
    ## min of synapsin and psd95
    img12 <- pmin(img1@.Data,img2@.Data)
    img12s <- img12/max(img12)
    return(list(egfp=img,egfp_max=f_max,egfp_mean=f_bar,
                egfp_max_mean=img0,egfp_max_mean_scl=img0s,
                synapsin=img1,synapsin_scl=img1s,psd=img2,psd_scl=img2s,
                synapse=img12,synapse_scl=img12s))
}

## plot using field::image -----------
my_heatmap <- function(dat,normalize=TRUE,truncate_thr=0.5) {
    if(normalize) {
        dat <- dat/max(dat)
    }
    dat[dat>truncate_thr] <- truncate_thr
    xx2 <- t(apply(dat, 1, rev))
    par(mar = c(.5,.5,.5,.5))
    fields::image.plot(xx2,axes=F,col = my_palette,legend.width=0.5)
    par(mar = c(2.5,2.5,2.5,2.5))
}

## remove out of boundary points ------------
my_within_boundary <- function(lst,Nx=512) {
    lst <- lst[lst[,1]>=1 & lst[,1]<=Nx &
                   lst[,2]>=1 & lst[,2]<=Nx,,drop=F]
}

## label a region ----------
my_label_region <- function(img,idx,idx2,idx3) {
    if(class(img)=='matrix') {
        img <- toRGB(as.Image(img))
    }
    if(img@colormode==0) {
        img <- toRGB((img))
    }
    if(dim(idx)[2]==dim(img)[2]) {
        idx <- which(idx>0,arr.ind = T)
    }
    img[cbind(idx[,1],idx[,2],1)] <- 1
    img[cbind(idx[,1],idx[,2],2)] <- 0
    img[cbind(idx[,1],idx[,2],3)] <- 0
    if(!missing(idx2)) {
        img[cbind(idx2[,1],idx2[,2],1)] <- 0
        img[cbind(idx2[,1],idx2[,2],2)] <- 1
        img[cbind(idx2[,1],idx2[,2],3)] <- 0
    }
    if(!missing(idx3)) {
        img[cbind(idx3[,1],idx3[,2],1)] <- 0
        img[cbind(idx3[,1],idx3[,2],2)] <- 0
        img[cbind(idx3[,1],idx3[,2],3)] <- 1
    }
    return(img)
}

my_label_region2 <- function(img,idx,idx2,idx3) {
    if(class(img)=='matrix') {
        img <- toRGB(as.Image(img))
    }
    if(img@colormode==0) {
        img <- toRGB((img))
    }
    if(dim(idx)[2]==dim(img)[2]) {
        idx <- which(idx>0,arr.ind = T)
    }
    img[cbind(idx[,1],idx[,2],1)] <- 0
    img[cbind(idx[,1],idx[,2],2)] <- 0
    img[cbind(idx[,1],idx[,2],3)] <- 1
    if(!missing(idx2)) {
        img[cbind(idx2[,1],idx2[,2],1)] <- 1
        img[cbind(idx2[,1],idx2[,2],2)] <- 0.5
        img[cbind(idx2[,1],idx2[,2],3)] <- 0
    }
    if(!missing(idx3)) {
        img[cbind(idx3[,1],idx3[,2],1)] <- 1
        img[cbind(idx3[,1],idx3[,2],2)] <- 0
        img[cbind(idx3[,1],idx3[,2],3)] <- 1
    }
    return(img)
}

my_label_region1 <- function(img,idx,idx2,idx3) {
    if(class(img)=='matrix') {
        img <- toRGB(as.Image(img))
    }
    if(img@colormode==0) {
        img <- toRGB((img))
    }
    if(dim(idx)[2]==dim(img)[2]) {
        idx <- which(idx>0,arr.ind = T)
    }
    img[cbind(idx[,1],idx[,2],1)] <- 1
    if(!missing(idx2)) {
        img[cbind(idx2[,1],idx2[,2],2)] <- 1
    }
    if(!missing(idx3)) {
        img[cbind(idx3[,1],idx3[,2],3)] <- 1
    }
    return(img)
}

## summation over surrounding pixels ------------
## !! must be square matrix
my_neib_sum <- function(dat) {
    Nx <- dim(dat)[1]
    Nx1 <- Nx-1
    dat_sum <- cbind(dat[,2:Nx],0) + cbind(0,dat[,1:Nx1]) +
        rbind(0,dat[1:Nx1,]) + rbind(dat[2:Nx,],0) +
        rbind(0,cbind(0,dat[1:Nx1,1:Nx1])) +
        rbind(0,cbind(dat[1:Nx1,2:Nx],0)) +
        rbind(cbind(dat[2:Nx,2:Nx],0),0) +
        rbind(cbind(0,dat[2:Nx,1:Nx1]),0)
    return(dat_sum)
}


