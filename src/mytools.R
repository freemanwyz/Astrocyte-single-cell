## label a region ----------
my_label_region <- function(idx,img) {
    if(class(img)=='matrix') {
        img <- toRGB(as.Image(img))
    }
    if(img1@colormode==0) {
        img <- toRGB((img))
    }
    if(dim(idx)[2]==dim(img)[2]) {
        idx <- which(idx>0,arr.ind = T)
    }
    img[cbind(idx[,1],idx[,2],1)] <- 1
    img[cbind(idx[,1],idx[,2],2)] <- 0
    img[cbind(idx[,1],idx[,2],3)] <- 0
    return(img)
}

## summation over surrounding pixels ------------
## !! must be square matrix
my_neib_sum <- function(dat) {
    Nx <- dim(dat)[1]
    Nx1 <- Nx-1
    dat_sum <- cbind(dat[,2:Nx],0) + cbind(0,dat[,1:Nx1]) +
        rbind(0,dat[1:Nx1,]) + rbind(dat[2:Nx,],0) +
        rbind(0,cbind(0,dat[1:Nx1,1:Nx1])) +
        rbind(0,cbind(dat[1:Nx1,2:Nx],0)) +
        rbind(cbind(dat[2:Nx,2:Nx],0),0) +
        rbind(cbind(0,dat[2:Nx,1:Nx1]),0)
    return(dat_sum)
}


