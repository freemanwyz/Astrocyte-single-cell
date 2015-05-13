## coculture syn1 spots ------------
source('./init_guilai_coculture.R')

img3 <- readImage(dat_Tuj1[[1]],type = 'tiff',all = T)
img3 <- channel(img3,'red')
# my_heatmap(img3,truncate_thr = 1)
img2 <- readImage(dat_Syn1[[1]],type = 'tiff',all = T)
img2 <- channel(img2,'red')
# my_heatmap(img2,truncate_thr = 1)
# display(img2)
# display(img3)

# img2s <- (img2>0.2)*1
# img3s <- (img3>0.2)*1
# display(img2s)
# display(img3s)

## process Tuj1 -----------
psize <- c(64,64)
npx <- dim(img3)[1]/psize[1]
npy <- dim(img3)[2]/psize[2]
# psize <- c(32,32)
# pcord <- c(6,7)
for(ii in seq(npx)) {
    for(jj in seq(npy)) {
        pcord <- c(ii,jj)
        rg1 <- seq((pcord[1]-1)*psize[1]+1,pcord[1]*psize[1])
        rg2 <- seq((pcord[2]-1)*psize[2]+1,pcord[2]*psize[2])
        img3_patch <- img3[rg1,rg2]
        fname <- paste0(ptmp,'/ctrl_neuron_ctrl_astro patch_',ii,'_',jj,'.tif')
        writeImage(img3_patch, files=fname, type='tiff', bits.per.sample=8)
    }
}
# my_heatmap(img3_patch,truncate_thr = 1)
