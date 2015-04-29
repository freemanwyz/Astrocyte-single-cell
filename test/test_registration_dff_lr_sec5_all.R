## relate synapse location to neural activity
## using AP and max-z

## init ----------
rm(list=ls())
source('./init.R')
activity_data <- 3

## neuron -----------------
if(activity_data==1) {  # neuron maxz
    img0 <- readImage(datfp$maxz_EGFP,type = 'tiff',all = T)
    img0 <- channel(img0,'green')
    ## heatmap
    xx0 <- img0
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
    show(sum(f_bar<0.0001))
    f_bar[f_bar<0.0001] <- 1
    rt0 <- (f_max-f_bar)/f_bar
    rt0[rt0<0] <- -1
    img0 <- as.Image(rt0/max(rt0))
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
    img0 <- rt0
    ## heatmap
    xx0 <- img0
    xx0 <- t(apply(xx0, 1, rev))
    par(mar = c(.5,.5,.5,.5))
    fields::image.plot(xx0,axes=F,col = my_palette,legend.width=0.5)
    par(mar = c(2.5,2.5,2.5,2.5))
}

## synapse --------------------
img1 <- readImage(datfp$synapsin[5],type = 'tiff',all = T)
img1 <- channel(img1,'gray')
img2 <- readImage(datfp$psd95[5],type = 'tiff',all = T)
img2 <- channel(img2,'gray')
img12 <- pmin(img1@.Data,img2@.Data)

## heatmap for img12 binarized
xx2 <- t(apply(img12, 1, rev))
par(mar = c(.5,.5,.5,.5))
fields::image.plot(xx2,axes=F,col = my_palette,legend.width=0.5)
par(mar = c(2.5,2.5,2.5,2.5))

## regression ----------------
idx <- which(img0>-1)
x <- img12[idx]
y <- img0[idx]
y_log <- log(y+0.01)
x_log <- log(x+0.01)
hist(x_log)
hist(y_log)

show(summary(lm(y~x)))
# show(summary(lm(y_log~x)))
show(summary(lm(y_log~x_log)))
par(mar = c(5,5,5,5))
plot(x_log,y_log,cex=0.1,xlab='log(x)',ylab='log(y)')



