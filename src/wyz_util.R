## preprocessing model fititng ----------
fn_prep <- function(x,dat,nk,nB,nC) {
    k <- x[nk]
    B <- x[nB]
    C <- x[nC]
    nPix <- dim(dat)[1]
    nTps <- dim(dat)[2]
    trg <- seq(nTps)
    f0 <- 1 + B*exp(-C*trg)
    ff <- 0
    for(t in trg) {
        ff <- ff + sum((dat[,t] - k*f0[t])^2)
    }
    return(ff)
}

gr_prep <- function(x,dat,nk,nB,nC) {
    k <- x[nk]
    B <- x[nB]
    C <- x[nC]
    nPix <- dim(dat)[1]
    nTps <- dim(dat)[2]
    trg <- seq(nTps)
    f0 <- 1 + B*exp(-C*trg)
    gr0 <- x*0
    ## k
    grk <- k*0
    for(t in trg) {
        grk <- grk + 2*(dat[,t] - k*f0[t])*(-1-B*exp(-C*t))
    }
    
    ## B
    grB <- 0
    for(t in trg) {
        grB <- grB + 2*sum( (dat[,t] - k*f0[t])*(-k*exp(-C*t)) )
    }
    
    ## C
    grC <- 0
    for(t in trg) {
        grC <- grC + 2*sum( (dat[,t] - k*f0[t])*(t*k*B*exp(-C*t)) )
    }
    # browser()
    return(c(grk,grB,grC))
}

gr_prep_numeric <- function(x,dat,nk,nB,nC) {
    gr0 <- x*0
    delta <- 1e-6
    for(ii in seq_along(x)) {
        xa <- x
        xb <- x
        xa[ii] <- xa[ii] + delta
        xb[ii] <- xb[ii] - delta
        gr0[ii] <- (fn_prep(xa,dat,nk,nB,nC) - fn_prep(xb,dat,nk,nB,nC))/delta/2
    }
    # browser()
    return(gr0)
}

## show many regions ------------------
show_many_regions <- function(bg, flst, cols_arr, pattern, minSize=50) {
    library(EBImage)
    d2 <- bg  # background image
    f2 <- flst  # ROI image list
    kk <- 1
    for(jj in seq_along(f2)) {
        d2x <- readImage(f2[jj])
        idx <- switch(pattern,
                      wyz=which(d2x[,,1]==1 & d2x[,,2]==0 & d2x[,,3]==1, arr.ind = T),
                      yinxue=which(d2x[,,1]>0 & d2x[,,2]==0 & d2x[,,3]==0, arr.ind = T)
        )
        # browser()
        if(dim(idx)[1]>minSize) {
            cat(jj,' ')
            d2[cbind(idx,1)] <- cols_arr[kk,1]
            d2[cbind(idx,2)] <- cols_arr[kk,2]
            d2[cbind(idx,3)] <- cols_arr[kk,3]
            kk <- kk + 1
            if(kk>dim(cols_arr)[1]) {
                kk <- 1
            }
        }
    }
    cat('\n')
    return(d2)
}

## generate region list from image files ------------------
gen_region_list <- function(flst, pattern, minSize=50) {
    library(EBImage)
    f2 <- flst  # ROI image list
    kk <- 1
    res <- list()
    res_name <- c()
    for(jj in seq_along(f2)) {
        d2x <- readImage(f2[jj])
        idx <- switch(pattern,
                      wyz=which(d2x[,,1]==1 & d2x[,,2]==0 & d2x[,,3]==1, arr.ind = T),
                      yinxue=which(d2x[,,1]>0 & d2x[,,2]==0 & d2x[,,3]==0, arr.ind = T)
        )
        if(dim(idx)[1]>minSize) {
            cat(jj,' ')
            res[[kk]] <- idx 
            res_name <- c(res_name,f2[jj])
            kk <- kk + 1
        }
    }
    cat('\n')
    names(res) <- res_name
    return(res)
}


## show a region ---------------
show_region <- function (corrMap, lstPixels, fname) {
    cmap1 <- corrMap
    cmap1[cmap1<0] <- 0
    img <- toRGB(as.Image(cmap1))
    idx <- do.call(rbind,lstPixels)
    idx <- cbind(idx,1)
    img[idx] <- 1
    if(missing(fname)) {
        display(transpose(img))
    } else {
        writeImage(transpose(img),files = fname, type = 'png')
    }
}

## downsample 2D or 3D ----------------
downsp_image <- function(dat,k,mthd='mean') {
    nx <- dim(dat)[1]
    ny <- dim(dat)[2]
    if(nx%%k>0 | ny%%k>0) {
        stop('k not supported yet\n')
    }
    nx1 <- nx/k
    ny1 <- ny/k
    res <- matrix(0, nx1, ny1)
    for(ii in seq(k)) {
        for(jj in seq(k)) {
            res <- res + dat[seq(ii,nx-k+ii,k),seq(jj,ny-k+jj,k)]
        }
    }
    res1 <- res/(k^2)
    return(res1)
}

downsp_video <- function(dat,k,mthd='mean') {
    nx <- dim(dat)[1]
    ny <- dim(dat)[2]
    if(nx%%k>0 | ny%%k>0) {
        stop('k not supported yet\n')
    }
    nx1 <- nx/k
    ny1 <- ny/k
    res <- array(0, c(nx1, ny1, dim(dat)[3]))
    for(ii in seq(k)) {
        for(jj in seq(k)) {
            res <- res + dat[seq(ii,nx-k+ii,k),seq(jj,ny-k+jj,k),]
        }
    }
    res1 <- res/(k^2)
    return(res1)
}


## fast correlation ----------------
fast_cor <- function(x,y) {
    nx <- length(x)
    ny <- length(y)
    ex <- sum(x)/nx
    ey <- sum(y)/ny
    ex2 <- sum(x^2)/nx
    e2x <- ex^2
    ey2 <- sum(y^2)/ny
    e2y <- ey^2
    sigx2 <- ex2 - e2x
    sigy2 <- ey2 - e2y
    if(sigx2>0 && sigy2>0){
        exy <- sum(x*y)/nx
        r <- (exy-ex*ey)/sqrt(sigx2*sigy2)
    }else{
        r <- 0
    }
    if(r>1) {
        r <- 1
        cat('rrr',r,'\n')
    }
    if(r< -1) {
        r <- -1
        cat('rrr',r,'\n')
    }
    return(r)
}

## find the pixels  ---------------
## in lstPixels but not in uniqUnit$lstUnitPixels
get_tmp_pix <- function(uniqUnit, lstPixels, nx=512, ny=512) {
    lstUnitPixels <- uniqUnit$lstUnitPixels
    idx1 <- do.call(rbind,lstUnitPixels)
    idx1s <- (idx1[,2]-1)*nx + idx1[,1]  # !! this function should know image size
    idx2 <- do.call(rbind,lstPixels)
    idx2s <- (idx2[,2]-1)*nx + idx2[,1]
    mat <- matrix(0,nx,nx)
    mat[idx2s] <- 1
    mat[idx1s] <- mat[idx1s] + 1
    res <- which(mat==1)
    return(res)
}

