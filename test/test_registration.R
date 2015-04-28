## relate synapse location to neural activity
## 40 stacks, single time point

rm(list=ls())
source('./init.R')
thr0 <- 0.03
thr1 <- 0.03
thr2 <- 0.03
plot_12 <- F
srg <- 1:10
ii <- 5
x <- c()
y <- c()
for(ii in srg) {
    cat(ii,' ')
    ## neuron
    img0 <- readImage(datfp$egfp[ii],type = 'tiff',all = T)
    img0 <- channel(img0,'gray')
    img0a <- img0
    img0a[img0<=thr0] <- 0
    img0a[img0>thr0] <- 1
    
    ## synapsin
    img1 <- readImage(datfp$synapsin[ii],type = 'tiff',all = T)
    img1 <- channel(img1,'gray')
    # display(img1)
    img1a <- img1
    img1a[img1<=thr1] <- 0
    img1a[img1>thr1] <- 1
    # display(img1a)
    
    ## psd95
    img2 <- readImage(datfp$psd95[ii],type = 'tiff',all = T)
    img2 <- channel(img2,'gray')
    # display(img2)
    img2a <- img2
    img2a[img2<=thr2] <- 0
    img2a[img2>thr2] <- 1
    # display(img2a)
    
    ## overlap
    img12a <- (img1a + img2a)/2
    if(plot_12) {
        idx1 <- which((img1a - img2a)>0,arr.ind = T)
        idx2 <- which((img2a - img1a)>0,arr.ind = T)
        idx12 <- which(img12a==1,arr.ind = T)
        img12 <- toRGB(img12a)
        img12[cbind(idx12[,1],idx12[,2],1)] <- 1
        img12[cbind(idx12[,1],idx12[,2],2)] <- 0
        img12[cbind(idx12[,1],idx12[,2],3)] <- 0
        img12[cbind(idx1[,1],idx1[,2],1)] <- 0
        img12[cbind(idx1[,1],idx1[,2],2)] <- 1
        img12[cbind(idx1[,1],idx1[,2],3)] <- 0
        img12[cbind(idx2[,1],idx2[,2],1)] <- 0
        img12[cbind(idx2[,1],idx2[,2],2)] <- 0
        img12[cbind(idx2[,1],idx2[,2],3)] <- 1
        display(img12)
    }
    
    ## regression for common parts
    img12b <- img12a
    img12b[img12b<1] <- 0
    idx_sel <- which(img0a>0)
    img12b[!idx_sel] <- 0
    x <- c(x,img12b[idx_sel])
    y <- c(y,img0[idx_sel])
}

idx0 <- x==0
idx1 <- x==1
y_log <- log(y)
# hist(y[idx0])
# hist(y[idx1])
hist(y_log[idx0])
hist(y_log[idx1])

summary(lm(y~x))
summary(lm(y_log~x))




