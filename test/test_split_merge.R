## on working

## init ----------
rm(list=ls())
source('./init.R')
thr_egfp_max_mean <- 0.03
hbound_region_grow <- 1.5

## load data
x <- my_load_dat()
img12 <- x$synapse
img12 <- img12/max(img12)
display(img12)

## find local maximum ------------
mat_val <- img12
ii <- 0
# while(min(mat_val)>0.05) {
while(ii<10) {
    seed <- which(mat_val==min(mat_val),arr.ind=T)[1,,drop=F]
    lst0 <- my_region_grow1(mat_val,seed,thr=hbound_region_grow)
    lst <- lst0$lst
    lst_all <- lst0$lst_all
    mat_val[lst_all] <- 100
    # tmp <- my_fill_hole(mat_val,lst_all)
    # mat_val <- tmp$mat
    display(my_label_region(mat_val,lst_all))
    ii <- ii + 1
}















