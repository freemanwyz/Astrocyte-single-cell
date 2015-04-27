## delta F/F for neuron activity
source('./init.R')

fname <- '40AP_40X'
fp <- paste0(pdat,'/',fname,'.tif')
img <- readImage(fp,type = 'tiff',all = T)
display(img)
img_neuron <- channel(img,'grey')
display(img_neuron)


