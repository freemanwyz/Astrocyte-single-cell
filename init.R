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

## parameters ---------
my_palette <- colorRampPalette(c('white',"yellow","green",'blue',"red"))(n = 299)

## find data folder from config file -------
if(!file.exists('~/wyz_dat')) {
    stop('Can not find configuration file\n')
}
tb <- read.table(file='~/wyz_dat',sep='=',stringsAsFactors = F,strip.white = T)
cfg <- as.list(tb[,2])
names(cfg) <- tb[,1] 





