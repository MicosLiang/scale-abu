area <- 1
dyarea <- 100
gyarea <- 4

'<-'(
  my_lm,
  function(x,y,d=1,draw=F)
  {
    ox <- x
    x <- sapply(0:d,function(v){return(x^v)})
    e <- try(w <- solve((t(x)%*%x))%*%(t(x)%*%y))
    if(typeof(e)=='character')
    {
      warning('not good dim')
      return(x,y,d-1,draw)
    }
    if(draw)
    {
      plot(ox,y,col=3,xlab='x',ylab='y',main=paste('lm: dim=',d))
      tx <- seq(min(ox)-10,max(ox)*1.2,0.1)
      x <- sapply(0:d,function(v){return(tx^v)})
      lines(tx, as.vector(x%*%w), col=2, type='l')
      legend('topright',legend=c('y=Σki*x^i',w))
    }
    return(w)
  }
)

'<-'(
  lm_t,
  function(x,y,w)
  {
    x <- cbind(1,x)
    yp <- x%*%w
    n <- dim(x)[1]
    p <- length(w)
    c <- solve(t(x)%*%x)
    tj <- numeric(p)
    s <- n-p-1
    t_sd <- sqrt((1/s)*sum((y-yp)^2))
    for(i in 1:p)
    {
      tj[i] <- -1 * w[i] / sqrt(c[i,i]*t_sd)
    }
    tj <- abs(tj)
    return(1-pt(tj,s))
  }
)

dat <- read.csv('./data/desert.csv')
dat[,3] <- tolower(dat[,3])
sps <- levels(factor(dat[,3]))
spv <- sapply(sps, function(sp){
  pdatid <- which(dat[,3]==sp)
  num <- length(pdatid)
  pdat <- dat[pdatid,]
  if(num!=1)
  {
    pdat <-  as.matrix(pdat[,4:6])
    v <- unlist(apply(pdat, 1, function(x){return(x[1] * x[2] * x[3])}))
  }
  else
  {
    pdat <- as.numeric(pdat[4:6])
    v <- pdat[1] * pdat[2] * pdat[3]
  }
  v <- mean(v)
  return(v)
})

dyang <- c('A','B','C','D','E','F')
d2yang <- c('M','N','O','P','Q','R','S','T')
gyang <- c('G','H','I','J','K','L','U','V')
ot <- c('OT')

ddat <- dat[which(dat[,1] %in% d2yang),]

iwant <- which(sps %in% levels(factor(ddat[,3])))

spnum <- sapply(sps, function(sp){
  pdatid <- which(ddat[,3]==sp)
  pdat <- ddat[pdatid,]
  if(length(pdatid)==1)
  {
    if(pdat[4]==0)
    {
      return(as.numeric(pdat[6]))
    }
  }
  else
  {
    if(all(pdat[,4]==0))
    {
      return(sum(pdat[,6]) / (dyarea * length(d2yang)) * area)
    }
  }
  return((length(pdatid) / (dyarea * length(d2yang))) * area)
})

spnum <- 100 * spnum / sum(spnum)

#ddat <- dat[which(dat[,1] %in% gyang),]

#spnum <- spnum + sapply(sps, function(sp){
#  pdatid <- which(ddat[,3]==sp)
#  pdat <- ddat[pdatid,]
#  if(length(pdatid)==1)
#  {
#    if(pdat[4]==0)
#    {
#      return(pdat[6]) 
#    }
#  }
#  else
#  {
#    if(all(pdat[,4]==0))
#    {
#      return(sum(pdat[,6]) / (gyarea * length(gyang)) * area)
#    }
#  }
#  return((length(pdatid) / (gyarea * length(gyang))) * area)
#})

spnum['cl'] <- 3
spnum['xzyhq'] <- 1
spnum['jje'] <- 17

spv <- ifelse(spv==0, 2, spv)
spv <- log(spv, 10)
spv <- spv[iwant]
spnum <- spnum[iwant]

w <- my_lm(spv, spnum, 3, T)
#lm_t(spv, spnum, w)
