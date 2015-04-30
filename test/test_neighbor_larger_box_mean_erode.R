## paired t-test
## relate synapse location to neural activity
## using AP and stack 5
## After thinning the mask, using larger box and take the mean

## init ----------
rm(list=ls())
source('./init.R')
neib_is_synapse <- T
erode <- F
thr0 <- 0.1
thr12 <- 0.12
thr12_low <- 0.06

## neuron -----------------
img <- readImage(datfp$ap,type = 'tiff',all = T)
img <- channel(img,'green')
f_max <- apply(img,c(1,2),max)
f_bar <- apply(img,c(1,2),mean)
rt0 <- f_max-f_bar
img0 <- rt0
display(img0)
img0a <- (img0>thr0)*1
display(img0a)

## thinning -----------------
im <- opening(closing(img0a))
im <- closing(img0a)
display(im)
im <- thinImage(im)
# im <- thinningIteration(im, 0)
# im <- thinningIteration(im, 1)
display(im)

# display(erode(img0a))
# display(closing(img0a))
display(erode(closing(img0a)))
if(erode) {
    img0a <- erode(closing(img0a))
}

## synapse --------------------
img1 <- readImage(datfp$maxz_synapsin,type = 'tiff',all = T)
img1 <- channel(img1,'gray')
img2 <- readImage(datfp$maxz_PSD95,type = 'tiff',all = T)
img2 <- channel(img2,'gray')
img12 <- pmin(img1@.Data,img2@.Data)
img12a <- (img12>thr12)*1
display(img12a)
img12a_low <- (img12>thr12_low)*1

## finding pairs ------------------
nb <- 5
nrg1 <- seq(-nb,nb)
nrg2 <- seq(-nb+1,nb-1)
ofst <- cbind(rbind(nrg1,0)+c(0,nb), rbind(nrg1,0)+c(0,-nb), 
              rbind(0,nrg2)+c(-nb,0), rbind(0,nrg2)+c(nb,0))
idx_sel <- which(img12a*img0a>0, arr.ind = T)
N_pix <- dim(idx_sel)[1]
mypairs_val <- list()
mypairs_pix <- list()
mypairs_neib <- list()
cat('N_pix:',N_pix,'\n')
jj <- 1
for(ii in seq(N_pix)) {
    if(ii%%100==0) cat(ii,' ')
    # cat(ii,' ')
    xx <- idx_sel[ii,]
    xx_neib <- t(ofst + xx)
    ## remove out of image neighbors
    xx_neib <- xx_neib[xx_neib[,1]>0 & xx_neib[,2]>0,,drop=F]
    xx_neib <- xx_neib[xx_neib[,1]<=512 & xx_neib[,2]<=512,,drop=F]
    if(length(xx_neib)==0) next
    ## remvoe pixels without activity
    xx_neib <- xx_neib[img0a[xx_neib]>0,,drop=F]
    if(length(xx_neib)==0) next
    ## remvoe pixels that are not synapse
    if(neib_is_synapse) {
        xx_neib <- xx_neib[img12a_low[xx_neib]>0,,drop=F]
        if(length(xx_neib)==0) next
    }
    ## select neighbors with lower value
    dif0 <- img12[xx[1],xx[2]] - img12[xx_neib]
    dif1 <- dif0[dif0>0]
    xx_neib <- xx_neib[dif0>0,,drop=F]
    if(length(xx_neib)==0) next
    ## gather results by mean
    mypairs_val[[jj]] <- c(img0[xx[1],xx[2]],mean(img0[xx_neib]))
    mypairs_pix[[jj]] <- as.vector(xx)
    mypairs_neib[[jj]] <- xx_neib
    jj <- jj+1
}

## analysis --------------
cat('N_pix_used:',length(mypairs_pix),'\n')
mypairs_val <- do.call(rbind,mypairs_val)
mypairs_pix <- do.call(rbind,mypairs_pix)
mypairs_neib <- do.call(rbind,mypairs_neib)

xx <- mypairs_pix[,1] + (mypairs_pix[,2]-1)*512
cat('Duplicated neighbors:')
show(table(table(xx)))

## show pairs
xx1 <- mypairs_pix
xx2 <- mypairs_neib
img0b <- toRGB(as.Image(img0/max(img0)))
img0b[cbind(xx1[,1],xx1[,2],1)] <- 1
img0b[cbind(xx1[,1],xx1[,2],2)] <- 0
img0b[cbind(xx1[,1],xx1[,2],3)] <- 0
img0b[cbind(xx2[,1],xx2[,2],1)] <- 0
img0b[cbind(xx2[,1],xx2[,2],2)] <- 1
img0b[cbind(xx2[,1],xx2[,2],3)] <- 0
display(img0b)

difx <- mypairs_val[,1]-mypairs_val[,2]
hist(difx,main='Pair difference distribution',xlab='pair difference')
show(t.test(mypairs_val[,1],mypairs_val[,2],paired = T))
cat('Positive ratio:',sum(difx>0)/length(difx),'\n')



