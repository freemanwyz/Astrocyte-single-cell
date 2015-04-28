## delta F/F for neuron activity
source('./init.R')

fname <- '40AP_40X'
fp <- paste0(pdat,'/',fname,'.tif')
img <- readImage(fp,type = 'tiff',all = T)
img <- channel(img,'green')
display(img)
N_tps <- dim(img)[3]

## first try -----------
f_bar <- apply(img,c(1,2),mean)
f_bar[f_bar==0] <- 2
dff <- img*0

for(ii in seq(N_tps)) {
    cat(ii,' ')
    dff[,,ii] <- (img[,,ii]-f_bar)/f_bar
}

mask_remove <- f_bar<0.01
dff_max <- apply(dff,c(1,2),max)
dff_max[mask_remove] <- 0
hmp <- dff_max/max(dff_max)
display(as.Image(hmp))

# f_bar_all <- mean(img)
# dff_all <- img/f_bar_all
# dff_all_max <- apply(dff_all,c(1,2),max)
# hmp_all <- dff_all_max/max(dff_all_max)
# display(as.Image(hmp_all))

plot(img[224,142,],type='l',col='blue',ylim=c(0,1))
lines(img[217,144,],col='red')

## bad --------------
f_max <- apply(img,c(1,2),max)
f_bar <- apply(img,c(1,2),mean)
f_bar[f_bar==0] <- 1
rt0 <- (f_max-f_bar)/f_bar
hmp1 <- rt0/max(rt0)
# hmp1[hmp1>0.6] <- 0.6
hmp1 <- hmp1/max(hmp1)
display(as.Image(hmp1))

## second try -----------
f_max <- apply(img,c(1,2),max)
f_bar <- apply(img,c(1,2),mean)
f_bar[f_bar<0.01] <- 1000
rt0 <- (f_max-f_bar)/f_bar
hmp1 <- rt0/max(rt0)
hmp1[hmp1>0.6] <- 0.6
hmp1 <- hmp1/max(hmp1)
display(as.Image(hmp1))

# idx <- which(hmp1>0.2,arr.ind = T)
# hmp2 <- toRGB(as.Image(hmp1))
# hmp2[cbind(idx[,1],idx[,2],1)] <- 1
# hmp2[cbind(idx[,1],idx[,2],2)] <- 0
# hmp2[cbind(idx[,1],idx[,2],3)] <- 0
# display(hmp2,title = '0.2')
# 
# idx <- which(hmp1>0.25,arr.ind = T)
# hmp2 <- toRGB(as.Image(hmp1))
# hmp2[cbind(idx[,1],idx[,2],1)] <- 1
# hmp2[cbind(idx[,1],idx[,2],2)] <- 0
# hmp2[cbind(idx[,1],idx[,2],3)] <- 0
# display(hmp2,title = '0.25')
# 
# idx <- which(hmp1>0.3,arr.ind = T)
# hmp2 <- toRGB(as.Image(hmp1))
# hmp2[cbind(idx[,1],idx[,2],1)] <- 1
# hmp2[cbind(idx[,1],idx[,2],2)] <- 0
# hmp2[cbind(idx[,1],idx[,2],3)] <- 0
# display(hmp2,title = '0.3')

idx <- which(hmp1>0.35,arr.ind = T)
hmp2 <- toRGB(as.Image(hmp1))
hmp2[cbind(idx[,1],idx[,2],1)] <- 1
hmp2[cbind(idx[,1],idx[,2],2)] <- 0
hmp2[cbind(idx[,1],idx[,2],3)] <- 0
display(hmp2,title = '0.35')

idx <- which(hmp1>0.4,arr.ind = T)
hmp2 <- toRGB(as.Image(hmp1))
hmp2[cbind(idx[,1],idx[,2],1)] <- 1
hmp2[cbind(idx[,1],idx[,2],2)] <- 0
hmp2[cbind(idx[,1],idx[,2],3)] <- 0
display(hmp2,title = '0.4')

idx <- which(hmp1>0.6,arr.ind = T)
hmp2 <- toRGB(as.Image(hmp1))
hmp2[cbind(idx[,1],idx[,2],1)] <- 1
hmp2[cbind(idx[,1],idx[,2],2)] <- 0
hmp2[cbind(idx[,1],idx[,2],3)] <- 0
display(hmp2,title = '0.5')


