## init -------
p0 <- list.files(pattern="*.R$", path="./src/", full.names=TRUE, ignore.case = T)
sapply(p0, FUN=source)

library(methods)
# library(tiff)
# library(png)
library(parallel)
library(gplots)
library(ggplot2)
# library(seriation)
# library(fields)
library(EBImage)
# library(xlsx)

## find data folder from config file -------
if(!file.exists('~/wyz_dat')) {
    stop('Can not find configuration file\n')
}

tb <- read.table(file='~/wyz_dat',sep='=',stringsAsFactors = F,strip.white = T)
cfg <- as.list(tb[,2])
names(cfg) <- tb[,1] 
pdat <- cfg$synapse_ucdavis_dat
ptmp <- cfg$synapse_registration_dump

## data files --------------
stack_idx <- 1:10
datfp <- list()
datfp$egfp <- paste0(pdat,'/tif/EGFP+synapsin+PSD 952.lsm - Ch1-T1 - C1 Z',stack_idx,' T1.tif')
datfp$synapsin <- paste0(pdat,'/tif/EGFP+synapsin+PSD 952.lsm - Ch2-T1 - C2 Z',stack_idx,' T1.tif')
datfp$psd95 <-  paste0(pdat,'/tif/EGFP+synapsin+PSD 952.lsm - Ch1-T2 - C3 Z',stack_idx,' T1.tif')
datfp$ap <- paste0(pdat,'/40AP_40X.tif')
datfp$maxz_EGFP <- paste0(pdat,'/MAX_EGFP.tif')
datfp$maxz_synapsin <- paste0(pdat,'/MAX_synapsin.tif')
datfp$maxz_PSD95 <- paste0(pdat,'/MAX_PSD952.tif')

## parameters ---------
my_palette <- colorRampPalette(c('white',"yellow","green",'blue',"red"))(n = 299)




