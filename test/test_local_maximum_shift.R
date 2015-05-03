## paired t-test
## relate synapse location to neural activity
## using shift

## init ----------
rm(list=ls())
source('./init.R')
thr_egfp_max_mean <- 0.03
hbound_region_grow <- 0.8
hbound_pair_sel <- 0.8
lbound <- 0.2
nn_max_ratio <- 0.1
val_thr <- 0.1

## load data
x <- my_load_dat()
xdat <- x$egfp_max_mean
img12 <- x$synapse * (xdat>thr_egfp_max_mean) * (x$egfp_mean<0.9) * 1
img12 <- img12/max(img12)
# my_heatmap(dat = img12,truncate_thr = 0.5)
# my_heatmap(dat = x$egfp_max_mean,truncate_thr = 1)

## setup
Nx <- dim(img12)[1]
rds <- 50
rds1 <- rds + 1
mat_ref <- matrix(0,rds+rds1,rds+rds1)
for(ii in seq(rds+rds1)) {
    for(jj in seq(rds+rds1)) {
        mat_ref[ii,jj] <- round(sqrt((ii-rds1)^2+(jj-rds1)^2))
    }
}

## find local maximum ------------
xegfp <- xdat
mat_val <- img12
mat_avail <- (mat_val>0)*1
res <- list()
pic <- img12
pic1 <- img12
ii <- 1
kk <- 0
while(max(mat_val)>val_thr) {
    kk <- kk + 1
    ## detect potential local maximum (ROI)
    seed <- which(mat_val==max(mat_val),arr.ind=T)[1,,drop=F]
    lst0 <- my_region_grow(mat_val,seed,thr=hbound_region_grow)
    res0 <- mean(mat_val[lst0])
    
    ## new centers
    # roi_center <- round(colMeans(lst0))
    circle_directions <- which(mat_ref==6,arr.ind = T) - rds1
    new_centers_mean <- rep(0,dim(circle_directions)[1])
    xx1_lst <- list()
    for(jj in seq(dim(circle_directions)[1])) {
        xx <- circle_directions[jj,]
        xx_proj <- colSums(t(lst0) * xx)
        idx_proj_max <- which.max(xx_proj)[1]
        idx_proj_min <- which.min(xx_proj)[1]
        dif0 <- lst0[idx_proj_max,] - lst0[idx_proj_min,]
        if(sum(dif0==0)==2) next
        xx1 <- t(t(lst0) + ceiling(dif0*1.2))
        xx1 <- my_within_boundary(xx1,Nx)
        xx1 <- xx1[mat_avail[xx1]>0,,drop=F]
        if(dim(xx1)[1]==0) next
        new_centers_mean[jj] <- mean(mat_val[xx1])
        xx1_lst[[jj]] <- xx1
    }
    
    if(sum(new_centers_mean)==0) {
        lst0a <- my_region_grow_n(dat = img12, lst = lst0, n = 2)
        mat_val[lst0a] <- 0
        mat_avail[lst0a] <- 0
        next
    }
    
    success <- 0
    nn <- 0
    nn_max <- ceiling(length(new_centers_mean)*nn_max_ratio)
    # nn_max <- 1
    while(max(new_centers_mean)>0 & nn<nn_max) {
        nn <- nn + 1
        ## opposite center to the max center
        idx1 <- which.max(new_centers_mean)[1]
        new_centers_mean[idx1] <- -1
        ## shifted lst0, remove small values
        xx1 <- xx1_lst[[idx1]]
        idx2a <- mat_val[xx1] >= lbound*mean(mat_val[lst0])
        idx2b <- mat_val[xx1] <= mean(mat_val[lst0])
        xx1 <- xx1[idx2a & idx2b,,drop=F]
        if(dim(xx1)[1]==0) next
        ## no overlap
        lst0s <- (lst0[,2]-1)*Nx+lst0[,1]
        xx1s <- (xx1[,2]-1)*Nx+xx1[,1]
        idxs <- xx1s %in% lst0s
        xx1 <- xx1[!idxs,,drop=F]
        if(dim(xx1)[1]==0) next
        ## criteria
        cr1 <- mean(mat_val[xx1]) < mean(mat_val[lst0])*hbound_pair_sel
        if( cr1 ) {  # target is smaller in synapse score
            res1 <- mean(xegfp[xx1])
            lst1 <- xx1
            success <- 1
            break
        }
    }
    
    ## update cache
    lst0a <- my_region_grow_n(dat = img12, lst = lst0, n = 2)
    mat_val[lst0a] <- 0
    mat_avail[lst0a] <- 0
    if(success) {
        if(res0<res1) {
            pic <- my_label_region2(pic,lst0,lst1)
            pic1 <- my_label_region2(pic1,lst0)
            pic1 <- my_label_region2(pic1,lst1)
        } else {
            pic <- my_label_region(pic,lst0,lst1)
            pic1 <- my_label_region(pic1,lst0)
            pic1 <- my_label_region(pic1,lst1)
        }
        lst1a <- my_region_grow_n(dat = img12, lst = lst1, n = 2)
        mat_val[lst1a] <- 0
        mat_avail[lst1a] <- 0
        res[[ii]] <- c(res0,res1)
        ii <- ii+1
        if(ii%%10==0) {
            cat(ii,' ')
        }
    }
}

display(pic)
# display(pic1)

## analysis ----------
res <- do.call(rbind,res)
res_dif <- res[,1] - res[,2]
hist(res[,1]/res[,2])
hist(res_dif)
show(t.test(res[,1],res[,2],paired = T))














