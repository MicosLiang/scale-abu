rm(list=ls())

'<-'(
  PCA,
  function(A, tag = NULL, center = T, normalize = F, d = 2, draw = T, show = T)
  {
    m <- dim(A)[1]
    n <- dim(A)[2]
    if(center)
    {
      A <- apply(A,2,function(x){return(x-mean(x))}) 
    }
    if(normalize)
    {
      A <- apply(A,2,function(x){return(scale(x))})
    }
    B <- 1/(m-1) * t(A) %*% A
    tmp <- eigen(B)
    values <- as.vector(tmp$values)
    disVaule <- values / sum(values)
    vectors <- as.matrix(tmp$vectors)
    vectors <- vectors[,1:d]
    C <- A %*% vectors
    if(d==2 && draw)
    {
      if(is.null(tag))
      {
        plot(C[,1],C[,2],xlab=paste('PC1(',disVaule[1],'%)'),ylab=paste('PC2(',disVaule[2],'%)'))
      }
      else
      {
        types <- levels(factor(tag))
        classes <- lapply(types, function(x){return(which(tag==x))})
        plot(C[classes[[1]],1],
             C[classes[[1]],2],
             col=2,
             xlab=paste('PC1(',round(disVaule[1]*100,2),'%)'),
             ylab=paste('PC2(',round(disVaule[2]*100,2),'%)'),
             xlim=c(min(C[,1]),max(C[,1])),
             ylim=c(min(C[,2]),max(C[,2])),
             pch=16)
        classes <- classes[-1]
        col_cnt <- 2
        for(each in classes)
        {
          col_cnt <- col_cnt + 1
          lines(C[each,1],C[each,2],col=col_cnt,pch=16,type='p')
        }
        if(show)
        {
          legend('topleft',legend = types, col=2:col_cnt, pch=16)  
        }
      }
    }
    return(C)
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

leaf_dat <- read.csv('./scale_abu_in_desert/leafs.csv')

envs <- levels(factor(leaf_dat[,1]))
jpeg('./scale_abu_in_desert/figure0.jpg', width=200*8, height=200*4)
par(mfcol=c(1,3))
lapply(envs, function(n){
  dat <- fildat(leaf_dat, 1, n)
  PCA(dat[,c(3:6,8:9)],draw = T,tag= dat[,2])
})

#lines(seq(0.1,2,0.1),sapply(seq(0.1,2,0.1),function(x){return(w[1]*x^w[2])}))