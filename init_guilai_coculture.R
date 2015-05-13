## init for guilai coculture -------
source('init.R')

pdat <- paste0(cfg$ucdavis_guilai_coculture,'/20150210_ctrl neuron_co_culture_tif/')
ptmp <- paste0(cfg$ucdavis_guilai_coculture,'/dump/')

## data files --------------
tb <- read.table(file='./test/coculture_flst.txt',sep='=',stringsAsFactors = F,strip.white = T)
cfg <- as.list(tb[,1])
dat_DAPI <- list()
dat_Syn1 <- list()
dat_Tuj1 <- list()
for(ii in seq_along(cfg)) {
    dat_DAPI[[ii]] <- paste0(pdat,'/',cfg[[ii]],' - Ch1-T1 - C1 Z1 T1.tif')
    dat_Syn1[[ii]] <- paste0(pdat,'/',cfg[[ii]],' - Ch2-T2 - C2 Z1 T1.tif')
    dat_Tuj1[[ii]] <- paste0(pdat,'/',cfg[[ii]],' - Ch3-T3 - C3 Z1 T1.tif')
}
names(dat_DAPI) <- as.vector(cfg)
names(dat_Syn1) <- as.vector(cfg)
names(dat_Tuj1) <- as.vector(cfg)
