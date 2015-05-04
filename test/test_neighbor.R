## paired t-test
## relate synapse location to neural activity
## using AP and stack 5

## init ----------
rm(list=ls())
source('./init.R')
activity_data <- 3
neib_is_synapse <- T
thr0 <- 0.03
# thr0 <- 0.1
# thr12 <- 0.05
thr12 <- 0.05
# thr12 <- 0.07
# thr12 <- 0.12
thr12_low <- 0.05
# thr12_low <- 0.07
stack_id <- 5

## neuron -----------------
if(activity_data==3) {  # neuron 40ap, f_max-f_mean
    img <- readImage(datfp$ap,type = 'tiff',all = T)
    img <- channel(img,'green')
    f_max <- apply(img,c(1,2),max)
    f_bar <- apply(img,c(1,2),mean)
    rt0 <- f_max-f_bar
    img0 <- rt0
    display(img0)
    img0a <- (img0>thr0)*1
    display(img0a)
}

## synapse --------------------
img1 <- readImage(datfp$maxz_synapsin,type = 'tiff',all = T)
# img1 <- readImage(datfp$synapsin[stack_id],type = 'tiff',all = T)
img1 <- channel(img1,'gray')
img2 <- readImage(datfp$maxz_PSD95,type = 'tiff',all = T)
# img2 <- readImage(datfp$psd95[stack_id],type = 'tiff',all = T)
img2 <- channel(img2,'gray')
img12 <- pmin(img1@.Data,img2@.Data)
img12a <- (img12>thr12)*1
display(img12a)
img12a_low <- (img12>thr12_low)*1
# display(img12a_low)

## finding pairs ------------------
nrg1 <- seq(-2,2)
nrg2 <- seq(-1,1)
ofst <- cbind(rbind(nrg1,0)+c(0,2), rbind(nrg1,0)+c(0,-2), 
                   rbind(0,nrg2)+c(-2,0), rbind(0,nrg2)+c(2,0))
idx_sel <- which(img12a*img0a>0,arr.ind = T)
N_pix <- dim(idx_sel)[1]
mypairs_val <- list()
# mypairs_syn <- list()
# mypairs_sin <- list()
# mypairs_psd <- list()
mypairs_pix <- list()
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
    ## sample a neighbor
    prob1 <- dif1/sum(dif1)
    idx_neib_sel <- sample(seq(length(dif1)),1,prob=prob1)
    neib_sel <- xx_neib[idx_neib_sel,]
    ## gather results
    mypairs_val[[jj]] <- c(img0[xx[1],xx[2]],img0[neib_sel[1],neib_sel[2]])
    # mypairs_syn[[jj]] <- c(img12[xx[1],xx[2]],img12[neib_sel[1],neib_sel[2]])
    # mypairs_sin[[jj]] <- c(img1[xx[1],xx[2]],img1[neib_sel[1],neib_sel[2]])
    # mypairs_psd[[jj]] <- c(img2[xx[1],xx[2]],img2[neib_sel[1],neib_sel[2]])
    mypairs_pix[[jj]] <- c(as.vector(xx),as.vector(neib_sel))
    jj <- jj+1
}
# cat('\n')

## analysis --------------
cat('N_pix_used:',length(mypairs_pix),'\n')
mypairs_val <- do.call(rbind,mypairs_val)
# mypairs_syn <- do.call(rbind,mypairs_syn)
# mypairs_sin <- do.call(rbind,mypairs_sin)
# mypairs_psd <- do.call(rbind,mypairs_psd)
mypairs_pix <- do.call(rbind,mypairs_pix)

xx <- mypairs_pix[,3] + (mypairs_pix[,4]-1)*512
cat('Duplicated neighbors:')
show(table(table(xx)))

## show pairs
xx1 <- mypairs_pix[,c(1,2)]
xx2 <- mypairs_pix[,c(3,4)]
img0b <- toRGB(as.Image(img0/max(img0)))
img0b[cbind(xx1[,1],xx1[,2],1)] <- 1
img0b[cbind(xx1[,1],xx1[,2],2)] <- 0
img0b[cbind(xx1[,1],xx1[,2],3)] <- 0
img0b[cbind(xx2[,1],xx2[,2],1)] <- 0
img0b[cbind(xx2[,1],xx2[,2],2)] <- 1
img0b[cbind(xx2[,1],xx2[,2],3)] <- 0
display(img0b)

difx <- mypairs_val[,1]-mypairs_val[,2]
# difsyn <- mypairs_syn[,1]-mypairs_syn[,2]
# difsin <- mypairs_sin[,1]-mypairs_sin[,2]
# difpsd <- mypairs_psd[,1]-mypairs_psd[,2]
# difmin <- mypairs_sin[,1]*mypairs_psd[,1]-mypairs_sin[,2]*mypairs_psd[,2]

hist(difx,main='Pair difference distribution',xlab='pair difference')
# hist(dify)
# plot(difx,dify)
show(t.test(mypairs_val[,1],mypairs_val[,2],paired = T))
cat('Positive ratio:',sum(difx>0)/length(difx),'\n')

# show(summary(lm(difx~difsyn)))
# show(summary(lm(difx~difsin)))
# show(summary(lm(difx~difpsd)))
# show(summary(lm(difx~difsin+difpsd)))
# show(summary(lm(difx~difsin+difpsd+difmin)))
# show(summary(lm(difx~difsin+difpsd+difsyn)))
# show(summary(lm(difx~difmin)))


