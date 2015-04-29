## paired t-test
## relate synapse location to neural activity
## using AP and stack 5

## init ----------
rm(list=ls())
source('./init.R')
activity_data <- 3
thr0 <- 0.1
thr12 <- 0.2

## neuron -----------------
if(activity_data==3) {  # neuron 40ap, f_max-f_mean
    img <- readImage(datfp$ap,type = 'tiff',all = T)
    img <- channel(img,'green')
    f_max <- apply(img,c(1,2),max)
    f_bar <- apply(img,c(1,2),mean)
    rt0 <- f_max-f_bar
    img0 <- rt0
    img0a <- (img0>thr0)*1
    display(img0a)
}

## synapse --------------------
img1 <- readImage(datfp$synapsin[5],type = 'tiff',all = T)
img1 <- channel(img1,'gray')
img2 <- readImage(datfp$psd95[5],type = 'tiff',all = T)
img2 <- channel(img2,'gray')
img12 <- pmin(img1@.Data,img2@.Data)
img12a <- (img12>thr12)*1
display(img12a)

## finding pairs ------------------
nrg1 <- seq(-2,2)
nrg2 <- seq(-1,1)
ofst <- cbind(rbind(nrg1,0)+c(0,2), rbind(nrg1,0)+c(0,-2), 
                   rbind(0,nrg2)+c(-2,0), rbind(0,nrg2)+c(2,0))
idx_sel <- which(img12a*img0a>0,arr.ind = T)
N_pix <- dim(idx_sel)[1]
mypairs_val <- list()
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
    mypairs_pix[[jj]] <- c(as.vector(xx),as.vector(neib_sel))
    jj <- jj+1
}

mypairs_val <- do.call(rbind,mypairs_val)
mypairs_pix <- do.call(rbind,mypairs_pix)

xx <- mypairs_pix[,3] + (mypairs_pix[,4]-1)*512
xx1 <- table(xx)

hist(mypairs_val[,1]-mypairs_val[,2])

t.test(mypairs_val[,1],mypairs_val[,2],paired = T)




