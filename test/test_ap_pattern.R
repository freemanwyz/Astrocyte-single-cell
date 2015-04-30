## activity pattern

rm(list=ls())
source('./init.R')
thr0 <- 0.1

img <- readImage(datfp$ap,type = 'tiff',all = T)
img <- channel(img,'green')
f_max <- apply(img,c(1,2),max)
f_bar <- apply(img,c(1,2),mean)
rt0 <- f_max-f_bar

idx0a <- which(rt0>thr0)
plot(f_bar[idx0a],f_max[idx0a],cex=0.1,xlab='mean',ylab='max')
plot(f_max[idx0a],f_bar[idx0a],cex=0.1,xlab='max',ylab='mean')
plot(f_bar[idx0a],rt0[idx0a],cex=0.1,xlab='mean',ylab='max-mean')

zz <- rt0[f_bar>0.05]
dff <- zz/f_bar[f_bar>0.05]
plot(f_bar[f_bar>0.05],dff,cex=0.1,xlab='mean',ylab='dff')
plot(rt0[f_bar>0.05],dff,cex=0.1,xlab='max-mean',ylab='dff')

img0c <- img
img0c[f_bar>0.5] <- 1
img0c[f_bar<=0.5] <- 0
display(img0c)
