rm(list=ls())

biodat <- read.csv('./scale_abu_in_desert/biomass.csv')
yangdat <- read.csv('./scale_abu_in_desert/desert.csv')

darea <- 100
garea <- 4

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