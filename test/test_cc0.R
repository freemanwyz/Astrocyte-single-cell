## coculture syn1 spots ---------
source('./init_guilai_coculture.R')

img3 <- readImage(dat_Tuj1[[1]],type = 'tiff',all = T)
img3 <- channel(img3,'red')
# my_heatmap(img3,truncate_thr = 1)
img2 <- readImage(dat_Syn1[[1]],type = 'tiff',all = T)
img2 <- channel(img2,'red')
# my_heatmap(img2,truncate_thr = 1)
display(img2)
display(img3)

img2s <- (img2>0.2)*1
img3s <- (img3>0.2)*1
display(img2s)
display(img3s)
