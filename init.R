## init -------
p0 <- list.files(pattern="*.R$", path="./src/", full.names=TRUE, ignore.case = T)
sapply(p0, FUN=source)

library(methods)
library(tiff)
library(png)
library(parallel)
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

