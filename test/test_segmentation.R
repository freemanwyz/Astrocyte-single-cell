## paired t-test
## relate synapse location to neural activity
## using AP and stack 5
## After thinning the mask, using larger box and take the mean

## init ----------
rm(list=ls())
source('./init.R')

## neuron -----------------
img <- readImage(datfp$ap,type = 'tiff',all = T)
img <- channel(img,'green')
dat <- img@.Data
f_max <- apply(dat,c(1,2),max)
f_bar <- apply(dat,c(1,2),mean)
img0 <- f_max-f_bar
display(img0>0.03)

nTps <- dim(dat)[3]
idx <- matrix(0,512,512)
idx[2:511,2:511] <- 1
idx <- which(idx>0)
dat1 <- dat
for(ii in seq(nTps)) {
    dat1[,,ii] <- dat1[,,ii] + matrix(1e-6*rnorm(512*512),512,512)
}
corrMap <- calCorrMap1(dat1,idx)
workZMap_tmp<-sqrt(nTps-3)*(log((1+corrMap)/(1-corrMap)))/2   
workZMap<-(workZMap_tmp-median(workZMap_tmp))/sqrt(0.516)

corrMap1 <- corrMap>0.2
corrMap1 <- my_neib_sum(corrMap1)>4 & corrMap1>0
display(corrMap1)

lstRegions<-genLstRegions(workZMap,30)


## synapse --------------------
img1 <- readImage(datfp$maxz_synapsin,type = 'tiff',all = T)
img1 <- channel(img1,'gray')
img1s <- img1/max(img1)
display(my_label_region(idx = img1s>150/256,img = img1s))

img2 <- readImage(datfp$maxz_PSD95,type = 'tiff',all = T)
img2 <- channel(img2,'gray')
img2s <- img2/max(img2)

img12 <- pmin(img1@.Data,img2@.Data)
img12 <- img12/max(img12)
display(my_label_region(idx = img1>1,img = img12))

img12e <- img12 * (img0>0.03)
xx2 <- t(apply(img12e, 1, rev))
xx2[xx2>0.5] <- 0.5
par(mar = c(.5,.5,.5,.5))
fields::image.plot(xx2,axes=F,col = my_palette,legend.width=0.5)
par(mar = c(2.5,2.5,2.5,2.5))

## super pixel
img12_4 <- downsp_image(img12,2)
display(img12>0.02)
display(img12_4>0.02)
img12c <- img12>0.01
img12d <- my_neib_sum(img12c)>2 & img12c>0
display(img12d)

## seeds -------------
dat <- img12@.Data
dat <- (dat>9/256)*1
dat <- my_neib_sum(dat)>3 & dat>0
display(my_label_region(idx = dat,img = img12))

dat <- img1@.Data
dat <- (dat>12/256)*1
dat <- my_neib_sum(dat)>3 & dat>0
display(my_label_region(idx = dat,img = img1))

X <- img1@.Data * dat

img1_bin <- (img1>7/256)*1
display(my_label_region(idx = img1_bin,img = img1))














