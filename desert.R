rm(list=ls())

'<-'(
  my_lm,
  function(x,y,d=1,draw=F,show=T,col=2)
  {
    ox <- x
    x <- sapply(0:d,function(v){return(x^v)})
    e <- try(w <- solve((t(x)%*%x))%*%(t(x)%*%y))
    if(typeof(e)=='character')
    {
      warning('not good dim')
      #return(my_lm(x,y,d-1,draw,show,col))
      return(F)
    }
    if(draw)
    {
      if(draw!=2)
      {
        plot(ox,y,col=3,xlab='x',ylab='y',main=paste('lm: dim=',d))  
      }
      tx <- seq(min(ox)-10,max(ox)*1.2,0.1)
      x <- sapply(0:d,function(v){return(tx^v)})
      lines(tx, as.vector(x%*%w), col=col, type='l')
      if(show){
        legend('topright',legend=c('y=Σki*x^i',w))
      }
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

'<-'(
  fildat,
  function(total_dat, col, ids)
  {
    if(length(ids > 1))
    {
      return(total_dat[which(total_dat[,col] %in% ids),])
    }
    return(total_dat[which(total_dat[,col]==ids),])
  }
)

biomass <- read.csv('./scale_abu_in_desert/biomass.csv')
biomass$v <- biomass$height * biomass$cw_max * biomass$cw_min * 0.001
sps <- levels(factor(biomass[,1]))
#plot(1,1,'n',xlim=c(min(biomass[,6]),max(biomass[,6])),ylim=c(min(biomass[,5]),max(biomass[,5])))
ws <- matrix(0,length(sps),2,dimnames=list(sps,c()))
plot(0,0,'n',xlim=c(-1,1),ylim=c(-1,1))
for(sp in 1:length(sps))
{
  spn <- sps[sp]
  pdatid <- which(biomass[,1]==spn)
  pdat <- biomass[pdatid,]
  if(all(pdat[,2]==0) || length(pdatid)==1)
  {
    next
  }
  if(all(pdat[,2]!=0) && all(pdat[,3]!=0))
  {
    #x <- scale(pdat[,2])
    #y <- scale(pdat[,5])
    x <- pdat[,7]
    y <- pdat[,6]
    w <- my_lm(x,y,1,F,F,sp)
  }
  else
  {
    #x <- scale(pdat[,6])
    #y <- scale(pdat[,5])
    x <- pdat[,2]
    y <- pdat[,6]
    w <- my_lm(x,y,1,F,F,sp)
    print(pdat)
  }
  #lines(x,y,'p',col=sp)
  ws[sp,] <- w
  #print(spn)
  #print(lm_t(x, y, w))
}
ws <- ws[-which(sps=='ltp'),]

area <- 1
dyarea <- 100
gyarea <- 4

dat <- read.csv('./scale_abu_in_desert/desert.csv')
dat[,3] <- tolower(dat[,3])
dat <- dat[-which(dat$sp_name=='cs' | dat$sp_name=='fs'),]
dat[which(dat$height!=0 & dat$cw_max == 0 & dat$cw_min == 0), 5:6] <- 1
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
  v <- median(v)
  return(v)
})

dyang <- c('A','B','C','D','E','F','W','X')
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
      return(sum(pdat[,6]) / (dyarea * length(dyang)) * area)
    }
  }
  return((length(pdatid) / (dyarea * length(dyang))) * area)
})

spnum <- spnum / sum(spnum)


spnum['cl'] <- 3
spnum['xzyhq'] <- 1
spnum['jje'] <- 17

spv <- ifelse(spv==0, 2, spv)
spv <- log(spv[iwant])
spnum <- spnum[iwant]

w <- my_lm(spv, spnum, 1, T)
lm_t(spv, spnum, w)

ddat <- dat[which(dat[,1] %in% gyang),]
iwant <- which(sps %in% levels(factor(ddat[,3])))
total_biomass <- unlist(apply(ddat, 1, function(x){
  if(x[3] %in% names(ws[,1]))
  {
    return(ws[x[3],1] + ws[x[3],2] * as.numeric(x[4]) * as.numeric(x[5]) * as.numeric(x[6]) * 0.001)
  }
  else
  {
    return(0)
  }
}))
total_biomass <- sapply(levels(factor(ddat[,3])), function(n){
  return(sum(total_biomass[which(ddat[,3]==n)]))
})
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
      return(sum(pdat[,6]) / (gyarea * length(gyang)) * area)
    }
  }
  return((length(pdatid) / (gyarea * length(gyang))) * area)
})
spnum <- spnum / sum(spnum)
spnum <- spnum[iwant]
spnum <- spnum[-which(total_biomass==0)]
total_biomass <- total_biomass[-which(total_biomass==0)]
total_biomass <- log(total_biomass)
lm_t(total_biomass, spnum, my_lm(total_biomass, spnum, 1, T))