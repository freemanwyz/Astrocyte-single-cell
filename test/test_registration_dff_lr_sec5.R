## relate synapse location to neural activity
## using AP and max-z

## init ----------
rm(list=ls())
source('./init.R')
thr0 <- 0.15
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
img1 <- readImage(datfp$synapsin[5],type = 'tiff',all = T)
img1 <- channel(img1,'gray')
img1a <- img1
img1a[img1<=thr1] <- 0
img1a[img1>thr1] <- 1
img2 <- readImage(datfp$psd95[5],type = 'tiff',all = T)
img2 <- channel(img2,'gray')
img2a <- img2
img2a[img2<=thr2] <- 0
img2a[img2>thr2] <- 1

## overlap
img12a <- img1a * img2a
img12 <- img1 * img2
img12 <- img12/max(img12)
img12[img0a==0] <- 0
idx_sel <- which(img0a>0)
x <- img12[idx_sel]
y <- img0[idx_sel]
x_bin <- img12a[idx_sel]

## heatmaps --------------
## heatmap for img12
# xx1 <- t(apply(img12, 1, rev))
# xx1[xx1>0.7] <- 0.7
# xx1[xx1<0.001] <- 0
# par(mar = c(.5,.5,.5,.5))
# fields::image.plot(xx1,axes=F,col = my_palette,legend.width=0.5)
# par(mar = c(2.5,2.5,2.5,2.5))

## heatmap for img12 binarized
x_bin1 <- img1a * img2a
xx2 <- t(apply(x_bin1, 1, rev))
par(mar = c(.5,.5,.5,.5))
fields::image.plot(xx2,axes=F,col = my_palette,legend.width=0.5)
par(mar = c(2.5,2.5,2.5,2.5))

## compare with xx0
xx2_bar <- 1-xx2
xx0_xx2 <- xx0 * xx2
xx0_xx2_bar <- xx0 * xx2_bar
par(mar = c(.5,.5,2.5,.5))
fields::image.plot(xx0_xx2,axes=F,col = my_palette,legend.width=0.5,main='pixels with synapse')
par(mar = c(2.5,2.5,2.5,2.5))
par(mar = c(.5,.5,2.5,.5))
fields::image.plot(xx0_xx2_bar,axes=F,col = my_palette,legend.width=0.5,main='pixels without synapse')
par(mar = c(2.5,2.5,2.5,2.5))

cat(sum(xx0_xx2>0),'\n')
cat(sum(xx0_xx2_bar>0),'\n')

## regression ----------------
y_log <- log(y)
x_log <- log(x+0.1)
# hist(x_log)
# hist(y_log)

hist(y_log[x_bin==0],main='Response with x=0')
hist(y_log[x_bin==1],main='Response with x=1')

# show(summary(lm(y~x)))
# show(summary(lm(y_log~x)))
# show(summary(lm(y_log~x_log)))
show(summary(lm(y_log~x_bin)))

# plot(x_bin,y_log,cex=0.1)



