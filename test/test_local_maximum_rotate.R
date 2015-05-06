## paired t-test
## relate synapse location to neural activity
## using AP and stack 5
## After thinning the mask, using larger box and take the mean

## init ----------
rm(list=ls())
source('./init.R')
thr_egfp_max_mean <- 0.03
# hbound <- 0.5
# hbound_region_grow <- 0.6
hbound_region_grow <- 0.75
hbound_pair_sel <- 0.9
# hbound <- 0.8
# lbound <- 0.5
# lbound <- 0.3
lbound <- 0.6
nn_max_ratio <- 0.1
val_thr <- 0.05

## load data
x <- my_load_dat()
# xdat <- x$egfp_max
xdat <- x$egfp_max_mean
img12 <- x$synapse * (xdat>thr_egfp_max_mean) * (x$egfp_mean<0.9) * 1
# img12 <- x$synapse * (x$egfp_max_mean>thr_egfp_max_mean)
img12 <- img12/max(img12)
# img12[img12>0.5] <- 0.5
# my_heatmap(dat = img12,truncate_thr = 0.5)
# my_heatmap(dat = x$egfp_max_mean,truncate_thr = 1)
# my_heatmap(dat = x$egfp_mean,truncate_thr = 0.3)

par(mar = c(2.5,2.5,2.5,2.5))
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
mat_val_center <- mat_val
mat_avail <- (mat_val>0)*1
# mat_avail <- ((xdat>thr_egfp_max_mean) & (mat_val>0))*1
res <- list()
pic <- img12
ii <- 1
kk <- 0
while(max(mat_val_center)>val_thr) {
    kk <- kk + 1
    ## detect potential local maximum (ROI)
    seed <- which(mat_val_center==max(mat_val_center),arr.ind=T)[1,,drop=F]
    lst0 <- my_region_grow(mat_val_center,seed,thr=hbound_region_grow)
    res0 <- mean(mat_val[lst0])
    
    ## roi property
    roi_center <- round(colMeans(lst0))
    dist <- sqrt(rowSums((t(t(lst0)-roi_center))^2))
    roi_radius <- round(max(dist)/1.5+1)
    # roi_radius <- ceiling(max(dist)/1.5)
    # gap <- round(mean(dist) + mean(dist)+1)
    gap <- round(max(dist) + mean(dist)+1)
    
    ## new centers
    idx <- which(mat_ref>=gap & mat_ref<=gap+1 ,arr.ind = T) - rds1
    new_centers <- t(t(idx) + roi_center)
    new_centers <- my_within_boundary(new_centers,Nx)
    idx <- mat_avail[new_centers]>0
    new_centers <- new_centers[idx,,drop=F]
    if(dim(new_centers)[1]==0) {
        lst0a <- my_region_grow_n(dat = img12, lst = lst0, n = 2)
        mat_val[lst0a] <- 0
        mat_avail[lst0a] <- 0
        mat_val_center[lst0a] <- 0
        next
    }
    
    ## search the best new center
    idx <- which(mat_ref<=roi_radius,arr.ind = T) - rds1
    new_centers_mean <- rep(0,dim(new_centers)[1])
    for(jj in seq(dim(new_centers)[1])) {
        xx <- new_centers[jj,]
        xx1 <- t(t(idx)+xx)
        xx1 <- my_within_boundary(xx1,Nx)
        xx1 <- xx1[mat_avail[xx1]>0,,drop=F]
        if(dim(xx1)[1]==0) next
        new_centers_mean[jj] <- mean(mat_val[xx1])
    }
    if(0) {
        idxc <- complex(real=idx[,1],imag=idx[,2])
        idx_arg <- Arg(idxc)
        plot(sort(idx_arg),new_centers_mean[order(idx_arg)],type='l')
        browser()
    }
    
    if( sum(colSums(abs(t(lst0)-c(237,150)))==0)>0) {
        browser()
    }
    
    ## use the opposite center that is neither too large nor too small
    success <- 0
    nn <- 0
    # nn_max <- ceiling(length(new_centers_mean)*nn_max_ratio)
    nn_max <- 1
    while(max(new_centers_mean)>0 & nn<nn_max) {
        nn <- nn + 1
        ## opposite center to the max center
        idx1 <- which.max(new_centers_mean)[1]
        new_centers_mean[idx1] <- -1
        big_center <- new_centers[idx1,,drop=F]
        # opp_center <- 2*roi_center - big_center
        opp_center <- big_center
        ## within the boundary, is available and is one of centers
        opp_center <- my_within_boundary(opp_center,Nx)
        if(dim(opp_center)[1]==0) next
        if(mat_avail[opp_center]==0) next
        idx2 <- which(new_centers[,1]==opp_center[1] & new_centers[,2]==opp_center[2])
        if(length(idx2)==0) next
        new_centers_mean[idx2] <- -1
        ## a circle around the opposite center
        xx1 <- t(t(idx)+as.vector(opp_center))
        xx1 <- my_within_boundary(xx1,Nx)
        xx1 <- xx1[mat_avail[xx1]>0,,drop=F]
        if(dim(xx1)[1]==0) next
        ## avoid overlap
        lst0s <- (lst0[,2]-1)*Nx+lst0[,1]
        xx1s <- (xx1[,2]-1)*Nx+xx1[,1]
        idxs <- xx1s %in% lst0s
        xx1 <- xx1[!idxs,,drop=F]
        if(dim(xx1)[1]==0) next
        ## remove small value pixels
        idx3 <- mat_val[xx1]>lbound*max(mat_val[xx1])
        # idx3 <- mat_val[xx1]>lbound*mean(mat_val[xx1])
        # idx3 <- mat_val[xx1]>lbound*mean(mat_val[lst0])
        # idx3 <- mat_val[xx1]<hbound*mat_val[seed] & mat_val[xx1]>lbound*mat_val[seed]
        xx1 <- xx1[idx3,,drop=F]
        if(dim(xx1)[1]==0) next
        cr1 <- mean(mat_val[xx1]) < mean(mat_val[lst0])*hbound_pair_sel
        # cr2 <- mean(mat_val[xx1]) > mean(mat_val[lst0])*lbound  # !! redundant
        # cr3 <- mean(mat_val[xx1]) > mean(mat_val[lst0])*(1+lbound)
        # if( cr1 & cr2 ) {  # target is smaller in synapse score
        if( cr1 ) {  # target is smaller in synapse score
            res1 <- mean(xegfp[xx1])
            lst1 <- xx1
            success <- 1
            break
        }
#         if( cr3 ) {  # target is larger
#             res1 <- res0
#             res0 <- mean(xegfp[xx1])
#             lst1 <- lst0
#             lst0 <- xx1
#             success <- 1
#             break
#         }
    }
    
    ## update cache
    lst0a <- my_region_grow_n(dat = img12, lst = lst0, n = 4)
    mat_val_center[lst0a] <- 0
    if(success) {
#         if(log(res0)-log(res1)>log(1.66)) {
#             pic <- my_label_region(pic,lst0)
#             pic <- my_label_region(pic,lst1)
#         } 
#         if(log(res0)-log(res1)<log(0.6)){
#             pic <- my_label_region2(pic,lst0)
#             pic <- my_label_region2(pic,lst1)
#         }
        
#         if(res0<res1-0.1) {
#             pic <- my_label_region2(pic,lst0)
#             pic <- my_label_region2(pic,lst1)
#         } 
#         if(res0>res1+0.1){
#             pic <- my_label_region(pic,lst0)
#             pic <- my_label_region(pic,lst1)
#         }
        
        if(res0<res1) {
            pic <- my_label_region2(pic,lst0,lst1)
#             pic <- my_label_region2(pic,lst0)
#             pic <- my_label_region2(pic,lst1)
        } else {
            pic <- my_label_region(pic,lst0,lst1)
#             pic <- my_label_region(pic,lst0)
#             pic <- my_label_region(pic,lst1)
        }
        # pic <- my_label_region(img12,lst0,lst1)
        # fname <- paste0(ptmp,'/rois/',ii,'.jpeg')
        # writeImage(pic, fname, quality=85)
        # display(pic)
        # display(my_label_region1(img12,lst0a,lst1a))
        lst1a <- my_region_grow_n(dat = img12, lst = lst1, n = 4)
        mat_val[lst0a] <- 0
        mat_val[lst1a] <- 0
        mat_avail[lst0a] <- 0
        mat_avail[lst1a] <- 0
        mat_val_center[lst1a] <- 0
        res[[ii]] <- c(res0,res1)
        ii <- ii+1
        if(ii%%10==0) {
            cat(ii,' ')
            # browser()
        }
    }
}

display(pic)

## analysis ----------
res <- do.call(rbind,res)
res_dif <- res[,1] - res[,2]
hist(res[,1]/res[,2])
par(mar = c(5,5,0.5,0.5))
hist(res_dif,main=NULL,xlab='x1-x2',breaks=50)
par(mar = c(2.5,2.5,2.5,2.5))
show(t.test(res[,1],res[,2],paired = T))














