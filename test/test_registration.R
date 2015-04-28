## relate synapse location to neural activity
source('./init.R')

srg <- 1:10
ii <- 5
for(ii in srg) {
    img1 <- readImage(datfp$synapsin[ii],type = 'tiff',all = T)
    img1 <- channel(img1,'gray')
    # display(img1)
    img2 <- readImage(datfp$psd95[ii],type = 'tiff',all = T)
    img2 <- channel(img2,'gray')
    # display(img2)
}