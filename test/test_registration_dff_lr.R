## relate synapse location to neural activity
## using AP and max-z
rm(list=ls())
source('./init.R')
thr0 <- 0.05
thr0a <- 0.02
thr1 <- 0.03
thr2 <- 0.03
activity_data <- 3

## neuron -----------------
if(activity_data==1) {  # neuron maxz
    img0 <- readImage(datfp$maxz_EGFP,type = 'tiff',all = T)
    img0 <- channel(img0,'green')
    img0a <- img0
    img0a[img0<=thr0] <- 0
    img0a[img0>thr0] <- 1
    display(img0a)
    ## heatmap
    xx0 <- img0
    xx0[xx0<thr0] <- 0
    xx0 <- t(apply(xx0, 1, rev))
    par(mar = c(.5,.5,.5,.5))
    fields::image.plot(xx0,axes=F,col = my_palette,legend.width=0.5)
    par(mar = c(2.5,2.5,2.5,2.5))
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
    ## heatmap
    xx0 <- rt0
    xx0[xx0>7] <- 7
    xx0 <- t(apply(xx0, 1, rev))
    par(mar = c(.5,.5,.5,.5))
    fields::image.plot(xx0,axes=F,col = my_palette,legend.width=0.5)
    par(mar = c(2.5,2.5,2.5,2.5))
}
if(activity_data==3) {  # neuron 40ap, f_max-f_mean
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
    ## heatmap
    xx0 <- img0
    xx0[xx0<thr0] <- 0
    xx0 <- t(apply(xx0, 1, rev))
    par(mar = c(.5,.5,.5,.5))
    fields::image.plot(xx0,axes=F,col = my_palette,legend.width=0.5)
    par(mar = c(2.5,2.5,2.5,2.5))
}

## synapse --------------------
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
idx_sel <- which(img0a>0)
img12 <- (img1/max(img1)+img2/max(img2))/2
img12[img0a==0] <- 0
## heatmap
xx1 <- t(apply(img12, 1, rev))
xx1[xx1>0.7] <- 0.7
par(mar = c(.5,.5,.5,.5))
fields::image.plot(xx1,axes=F,col = my_palette,legend.width=0.5)
par(mar = c(2.5,2.5,2.5,2.5))

x <- img12[idx_sel]
y <- img0[idx_sel]

## regression ----------------
y_log <- log(y)
x_log <- log(x+0.1)
hist(x_log)
hist(y_log)

show(summary(lm(y~x)))
show(summary(lm(y_log~x)))
show(summary(lm(y_log~x_log)))

plot(x_log,y_log,cex=0.1)



