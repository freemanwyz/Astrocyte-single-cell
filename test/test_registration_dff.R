## relate synapse location to neural activity
## using AP and max-z
rm(list=ls())
source('./init.R')
thr0 <- 0.13
thr0a <- 0.02
thr1 <- 0.03
thr2 <- 0.03
activity_data <- 2
use_cor <- F

## charX from 40AP
if(use_cor) {
    img <- readImage(datfp$ap,type = 'tiff',all = T)
    img <- channel(img,'green')
    N_tps <- dim(img)[3]
    f_bar <- apply(img,c(1,2),mean)
    idx <- which(f_bar>0.05,arr.ind = T)
    N_pix <- dim(idx)[1]
    dat <- matrix(0,N_pix,N_tps)
    for(ii in seq(N_tps)) {
        cat(ii,' ')
        dat[,ii] <- img[cbind(idx[,1],idx[,2],ii)]
    }
    charX <- colMeans(dat)
    
    dat0 <- matrix(0,512*512,N_tps)
    idx <- which(f_bar>=0,arr.ind = T)
    for(ii in seq(N_tps)) {
        cat(ii,' ')
        dat0[,ii] <- img[cbind(idx[,1],idx[,2],ii)] + 1e-6*rnorm(512*512)
    }
    res_cor <- cor(charX,t(dat0))
    cor_map <- matrix(0,512,512)
    idx_cor <- idx[res_cor>0.9,]
    cor_map[idx_cor] <- 1
    display(as.Image(cor_map))
}

## neuron
if(activity_data==1) {  # neuron maxz
    img0 <- readImage(datfp$maxz_EGFP,type = 'tiff',all = T)
    img0 <- channel(img0,'green')
    img0a <- img0
    img0a[img0<=thr0] <- 0
    img0a[img0>thr0] <- 1
    display(img0a)
}
if(activity_data==2) {  # neuron 40ap, max(dff)
    img <- readImage(datfp$ap,type = 'tiff',all = T)
    img <- channel(img,'green')
    f_max <- apply(img,c(1,2),max)
    f_bar0 <- apply(img,c(1,2),mean)
    f_bar <- f_bar0
    f_bar[f_bar<thr0a] <- 1
    rt0 <- (f_max-f_bar)/f_bar
    rt0[rt0<0] <- 0
    img0 <- as.Image(rt0/max(rt0))
    img0a <- img0
    img0a[img0<=thr0] <- 0
    img0a[img0>thr0] <- 1
    display(img0a)
    xx0 <- rt0
    xx0[xx0>7] <- 7
    xx0 <- t(apply(xx0, 1, rev))
    my_palette <- colorRampPalette(c('white',"yellow","green",'blue',"red"))(n = 299)
    fields::image.plot(xx0,axes=F,col = my_palette,legend.width=0.5)
    
}
if(activity_data==3) {  # neuron 40ap, max
    img <- readImage(datfp$ap,type = 'tiff',all = T)
    img <- channel(img,'green')
    f_max <- apply(img,c(1,2),max)
    img0 <- as.Image(f_max)
    img0a <- img0
    img0a[img0<=thr0] <- 0
    img0a[img0>thr0] <- 1
    display(img0a)
}
if(activity_data==4) {  # neuron 40ap, max/min
    img <- readImage(datfp$ap,type = 'tiff',all = T)
    img <- channel(img,'green')
    f_max <- apply(img,c(1,2),max)
    f_min <- apply(img,c(1,2),min)
    f_min[f_min==0] <- 1
    # rt0 <- (f_max-f_min)/f_min
    rt0 <- f_max/f_min
    img0 <- as.Image(rt0/max(rt0))
    img0a <- img0
    img0a[img0<=thr0] <- 0
    img0a[img0>thr0] <- 1
    display(img0a)
}
if(activity_data==5) {  # neuron 40ap, f_max-f_mean
    img <- readImage(datfp$ap,type = 'tiff',all = T)
    img <- channel(img,'green')
    f_max <- apply(img,c(1,2),max)
    f_bar <- apply(img,c(1,2),mean)
    rt0 <- f_max-f_bar
    rt0[rt0<0] <- 0
    img0 <- as.Image(rt0/max(rt0))
    img0a <- img0
    img0a[img0<=thr0] <- 0
    img0a[img0>thr0] <- 1
    display(img0a)
}

## synapsin
img1 <- readImage(datfp$maxz_synapsin,type = 'tiff',all = T)
img1 <- channel(img1,'gray')
img1a <- img1
img1a[img1<=thr1] <- 0
img1a[img1>thr1] <- 1

## psd95
img2 <- readImage(datfp$maxz_PSD952,type = 'tiff',all = T)
img2 <- channel(img2,'gray')
img2a <- img2
img2a[img2<=thr2] <- 0
img2a[img2>thr2] <- 1

## overlap
img12a <- (img1a + img2a)/2
img12b <- img12a
img12b[img12b<1] <- 0
if(use_cor) {
    idx_sel <- which(img0a>0 & cor_map>0)
} else {
    idx_sel <- which(img0a>0)
}
img12b[!idx_sel] <- 0
x <- img12b[idx_sel]
y <- img0[idx_sel]

## regression for common parts
idx0 <- x==0
idx1 <- x==1
y_log <- log(y)
hist(y_log[idx0])
hist(y_log[idx1])

show(summary(lm(y~x)))
show(summary(lm(y_log~x)))




