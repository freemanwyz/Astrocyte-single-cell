## paired test or not
rm(list=ls())
source('init.R')
cond <- 1  # 1: no within group difference
balanced <- F
par(mar = c(5,5,2,2))

## generate data --------------
set.seed(889)
dat <- as.data.frame(matrix(0,150,3))
colnames(dat) <- c('synapse','calcium','group')
dat[1:50,1] <- 120 + 0.5*rnorm(50)
dat[51:100,1] <- 80 + 0.5*rnorm(50)
dat[101:125,1] <- 100 + 0.3*rnorm(25)
dat[126:150,1] <- 60 + 0.3*rnorm(25)
if(cond==1) {
    dat[1:100,2] <- 20 + rnorm(100)
    dat[101:150,2] <- 15 + 0.5*rnorm(50)
} else {
    dat[1:50,2] <- 19.75 + rnorm(50)
    dat[51:100,2] <- 20.25 + rnorm(50)
    dat[101:125,2] <- 14.75 + 0.5*rnorm(25)
    dat[126:150,2] <- 15.25 + 0.5*rnorm(25)
}
# dat[c(1:50,101:125),3] <- 1
# dat[c(51:100,126:150),3] <- 0
dat[c(1:50,101:125),3] <- '1'
dat[c(51:100,126:150),3] <- '0'

## two group tests -----------
if(!balanced) {
    dat <- dat[dat[,1]>70,]
    x1a <- dat[c(1:50,101:125),2]
    x2a <- dat[c(51:100),2]
    show(t.test(x1a,x2a))
    show(wilcox.test(x1a,x2a))
    x1b <- dat[c(1:50,101:110),2]
    x2b <- dat[c(51:100,111:120),2]
    show(t.test(x1b,x2b,paired = T))
    show(wilcox.test(x1b,x2b,paired = T))
    x12 <- data.frame(x1b-x2b)
} else {
    x1 <- dat[c(1:50,101:125),2]
    x2 <- dat[c(51:100,126:150),2]
    show(t.test(x1,x2))
    show(wilcox.test(x1,x2))
    show(t.test(x1,x2,paired = T))
    show(wilcox.test(x1,x2,paired = T))
    x12 <- data.frame(x1-x2)
}
colnames(x12) <- 'diff0'
pic1 <- ggplot(dat, aes(x = calcium, fill = group)) + 
    geom_histogram(binwidth=.5, alpha=.5, position="identity") + 
    xlab('Calcium score')
pic2 <- ggplot(x12, aes(x = diff0)) + 
    geom_histogram(binwidth=.5, alpha=.5, position="identity") + 
    xlab('Calcium score difference')

## linear regression --------------
show(summary(lm(dat[,2]~dat[,1])))
plot(dat[,1],dat[,2],xlab='Synapse score',ylab='Calcium score')
abline(lm(dat[,2]~dat[,1]),col='red')
show(summary(lm(dat[,2]~(dat[,3]==0))))
plot(dat[,3]==1,dat[,2],xlab='group',ylab='Calcium score')
abline(lm(dat[,2]~(dat[,3]==1)),col='red')






