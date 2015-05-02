## paired t-test
## relate synapse location to neural activity
## using AP and stack 5
## After thinning the mask, using larger box and take the mean

## init ----------
rm(list=ls())
source('./init.R')
thr_egfp_max_mean <- 0.03

## load -----------------
x <- my_load_dat()
img12 <- x$synapse * (x$egfp_max_mean>thr_egfp_max_mean)
img12 <- img12/max(img12)
img12[img12>0.5] <- 0.5
my_heatmap(dat = img12,truncate_thr = 1)

## find local maximum ------------
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

## detect potential local maximum (ROI)
seed <- which(img12==max(img12),arr.ind=T)[1,]
lst <- my_region_grow(img12,seed)

## roi property
roi_center <- round(colMeans(lst))
dist <- sqrt(rowSums((t(t(lst)-roi_center))^2))
roi_radius <- ceiling(max(dist)/1.5)
gap <- round(max(dist) + mean(dist))

## new centers
idx <- which(mat_ref==gap,arr.ind = T) - rds1
new_centers <- t(t(idx) + roi_center)
new_centers <- my_within_boundary(new_centers,Nx)
idx <- x$egfp_max_mean_scl[new_centers]>thr_egfp_max_mean
new_centers <- new_centers[idx,,drop=F]

## search the best new center
idx <- which(mat_ref<=roi_radius,arr.ind = T) - rds1
for(ii in seq(dim(new_centers)[1])) {
    
}

display(my_label_region(img12,lst,new_centers))














