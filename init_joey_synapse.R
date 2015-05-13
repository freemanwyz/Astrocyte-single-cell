source('init.R')

pdat <- paste0(cfg$synapse_ucdavis,'/data_joey_synapse/')
# ptmp <- cfg$synapse_registration_dump

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