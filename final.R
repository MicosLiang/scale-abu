rm(list=ls())
#library(ggplot2)
#library(tableone)
#library(deSolve)

biodat <- read.csv('./scale_abu_in_desert/biomass.csv')
yangdat <- read.csv('./scale_abu_in_desert/desert.csv')
yangdat[,3] <- tolower(yangdat[,3])
yangdat[which(yangdat$height!=0 & yangdat$cw_max == 0 & yangdat$cw_min == 0), 5:6] <- 1
yangdat <- yangdat[-which(yangdat$sp_name=='cs' | yangdat$sp_name=='fs'),]

darea <- 100
garea <- 4
yangnum <- 8

dyang <- c('A','B','C','D','E','F','W','X')
dyang2 <- c('M','N','O','P','Q','R','S','T')
gyang <- c('G','H','I','J','K','L','U','V')

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

'<-'(
  calS,
  function(x)
  {
    x <- as.numeric(x[4:6])
    return(x[2] * x[3] / 4 * pi)
  }
)

'<-'(
  calV,
  function(x)
  {
    x <- as.numeric(x[4:6])
    return(x[1] * x[2] * x[3] / 4 * pi)
  }
)

'<-'(
  juan, #一个类似卷积的操作
  function(dm,func)
  {
    len <- dim(dm)[1]
    ans <- NULL
    ans_scale <- character()
    ans_id <- character()
    al_cnt <- 1
    for(scale in 1:(len-1))
    {
      square <- (scale+1)^2
      end <- len - scale
      cnt <- 1
      for(posx in end:1)
      {
        for(posy in 1:end)
        {
          edgex <- posx + scale
          edgey <- posy + scale
          ans[al_cnt] <- func(dm[posx:edgex,posy:edgey])
          #ans_scale[al_cnt] <- paste('scale',scale,sep="",collaspe="")
          ans_scale[al_cnt] <- square
          #ans_id[al_cnt] <- paste('i',cnt,sep="",collapse = "")
          cnt <- cnt + 1
          al_cnt <- al_cnt + 1
        }
      }
    }
    ans <- cbind(as.numeric(ans_scale),as.numeric(ans))
    return(ans)
  }
)

#分析群落结构

#逻辑斯蒂精度
qun_step <- 0.2

#采样地1
dat1 <- fildat(yangdat, 1, dyang)
#所含物种
sps1 <- levels(factor(dat1[,3]))

#种群特征
jpeg('./scale_abu_in_desert/figure1.jpg', width=200*8, height=200*4)
#split.screen(c(ceiling(length(sps1)^0.5), ceiling(length(sps1)^0.5)))
par(mfcol=c(2,3))
sapply(1:length(sps1), function(n){
  spn <- sps1[n]
  spdat <- fildat(dat1, 3, spn)
  spv <- log(as.numeric(apply(spdat, 1, calV)))
  #screen(n)
  m <- seq(min(spv),max(spv),qun_step)
  n <- sapply(seq(min(spv),max(spv),qun_step), function(x){
    return(length(which(spv<=x)))
  })
  plot(seq(min(spv),max(spv),qun_step), sapply(seq(min(spv),max(spv),qun_step), function(x){
    return(length(which(spv<=x)))
  }), main = spn, xlab = 'log(V)', ylab = 'ind num', col = 4)
  '<-'(
    tmpfunc,
    function()
    {
      df <- as.data.frame(cbind(m,n))
      SS <- getInitial(n ~ SSlogis(m, alpha, xmid, scale), data = df)
      K_start <- SS["alpha"]
      R_start <- 1/SS["scale"]
      N0_start <- SS["alpha"]/(exp(SS["xmid"]/SS["scale"])+1)
      
      log_formula <- formula(n ~ K*N0*exp(R*m)/(K + N0*(exp(R*m) - 1)))
      formu<-nls(log_formula, start = list(K = K_start, R = R_start, N0 = N0_start))
      K <- coef(formu)[1]
      R <- coef(formu)[2]
      N0 <- coef(formu)[3]
      lines(min(m):max(m),sapply(min(m):max(m),function(x){
        return(K*N0*exp(R*x)/(K + N0*(exp(R*x) - 1)))
      }), col = 3, lty = 2)
    }
  )
  #try(tmpfunc())
})
dev.off()

#采样地2
dat2 <- fildat(yangdat, 1, dyang2)
#所含物种
sps2 <- levels(factor(dat2[,3]))

#种群特征
jpeg('./scale_abu_in_desert/figure2.jpg', width=200*ceiling(length(sps2)^0.5), height=200*ceiling(length(sps2)^0.5))
par(mfcol=c(ceiling(length(sps1)^0.5),ceiling(length(sps1)^0.5)))
sapply(1:length(sps2), function(n){
  spn <- sps2[n]
  spdat <- fildat(dat2, 3, spn)
  spv <- log(as.numeric(apply(spdat, 1, calV)))
  m <- seq(min(spv),max(spv),qun_step)
  n <- sapply(seq(min(spv),max(spv),qun_step), function(x){
    return(length(which(spv<=x)))
  })
  plot(seq(min(spv),max(spv),qun_step), sapply(seq(min(spv),max(spv),qun_step), function(x){
    return(length(which(spv<=x)))
  }), main = spn, xlab = 'log(V)', ylab = 'ind num', col = 7)
  '<-'(
    tmpfunc,
    function()
    {
      df <- as.data.frame(cbind(m,n))
      SS <- getInitial(n ~ SSlogis(m, alpha, xmid, scale), data = df)
      K_start <- SS["alpha"]
      R_start <- 1/SS["scale"]
      N0_start <- SS["alpha"]/(exp(SS["xmid"]/SS["scale"])+1)
      
      log_formula <- formula(n ~ K*N0*exp(R*m)/(K + N0*(exp(R*m) - 1)))
      formu<-nls(log_formula, start = list(K = K_start, R = R_start, N0 = N0_start))
      K <- coef(formu)[1]
      R <- coef(formu)[2]
      N0 <- coef(formu)[3]
      lines(min(m):max(m),sapply(min(m):max(m),function(x){
        return(K*N0*exp(R*x)/(K + N0*(exp(R*x) - 1)))
      }), col = 1, lty = 2)
    }
  )
  #try(tmpfunc())
})
dev.off()

#采样地3
dat3 <- fildat(yangdat, 1, gyang)
#所含物种
sps3 <- levels(factor(dat3[,3]))

#种群特征
#jpeg('./scale_abu_in_desert/figure3.jpg', width=200*ceiling(length(sps3)^0.5), height=200*ceiling(length(sps3)^0.5))
#par(mfcol=c(ceiling(length(sps3)^0.5),ceiling(length(sps3)^0.5)))
jpeg('./scale_abu_in_desert/figure3.jpg', width=200*8, height=200*8)
par(mfcol=c(4,4))
sapply(1:length(sps3), function(n){
  spn <- sps3[n]
  spdat <- fildat(dat3, 3, spn)
  spv <- log(as.numeric(apply(spdat, 1, calV)))
  m <- seq(min(spv),max(spv),qun_step)
  n <- sapply(seq(min(spv),max(spv),qun_step), function(x){
    return(length(which(spv<=x)))
  })
  plot(m, n, main = spn, xlab = 'log(V)', ylab = 'ind num', col = 6)
  '<-'(
    tmpfunc,
    function()
    {
      df <- as.data.frame(cbind(m,n))
      SS <- getInitial(n ~ SSlogis(m, alpha, xmid, scale), data = df)
      K_start <- SS["alpha"]
      R_start <- 1/SS["scale"]
      N0_start <- SS["alpha"]/(exp(SS["xmid"]/SS["scale"])+1)
      
      log_formula <- formula(n ~ K*N0*exp(R*m)/(K + N0*(exp(R*m) - 1)))
      formu<-nls(log_formula, start = list(K = K_start, R = R_start, N0 = N0_start))
      K <- coef(formu)[1]
      R <- coef(formu)[2]
      N0 <- coef(formu)[3]
      lines(min(m):max(m),sapply(min(m):max(m),function(x){
        return(K*N0*exp(R*x)/(K + N0*(exp(R*x) - 1)))
      }), col = 5, lty = 2)
    }
  )
  #try(tmpfunc())
})
dev.off()


#优势度

#盖度
space1 <- sapply(1:length(sps1), function(n){
  spn <- sps1[n]
  spdat <- fildat(dat1, 3, spn)
  spspace <- (sum(as.numeric(apply(spdat, 1, calS))) / 10000) / (darea * yangnum)
  return(spspace)
})
names(space1) <- sps1

space2 <- sapply(1:length(sps2), function(n){
  spn <- sps2[n]
  spdat <- fildat(dat2, 3, spn)
  spspace <- (sum(as.numeric(apply(spdat, 1, calS))) / 10000) / (darea * yangnum)
  return(spspace)
})
names(space2) <- sps2

space3 <- sapply(1:length(sps3), function(n){
  spn <- sps3[n]
  spdat <- fildat(dat3, 3, spn)
  spspace <- (sum(as.numeric(apply(spdat, 1, calS))) / 10000) / (garea * yangnum)
  return(spspace)
})
names(space3) <- sps3

#shannon-winner
spnum1 <- sapply(1:length(sps1), function(n){
  spn <- sps1[n]
  return(length(which(dat1[,3]==spn)))
})
names(spnum1) <- sps1
pis1 <- spnum1/sum(spnum1)
shannon1 <- -sum(pis1 * log(pis1))

spnum2 <- sapply(1:length(sps2), function(n){
  spn <- sps2[n]
  return(length(which(dat2[,3]==spn)))
})
names(spnum2) <- sps2
pis2 <- spnum2/sum(spnum2)
shannon2 <- -sum(pis2 * log(pis2))

spnum3 <- sapply(1:length(sps3), function(n){
  spn <- sps3[n]
  return(length(which(dat3[,3]==spn)))
})
names(spnum3) <- sps3
pis3 <- spnum3/sum(spnum3)
shannon3 <- -sum(pis3 * log(pis3))

#t.test(c(shannon1,shannon2,shannon3))

#simpson
simpson1 <- 1 - sum(pis1^2)

simpson2 <- 1 - sum(pis2^2)

simpson3 <- 1 - sum(pis3^2)

#t.test(c(simpson1,simpson2,simpson3))

#多样性分解
alpha_div1 <- length(sps1)
alpha_div2 <- length(sps2)
alpha_div3 <- length(sps3)

gamma_div <- length(union(union(sps1, sps2), sps3))

#beta_div <- 1 - mean(c(alpha_div1,alpha_div2,alpha_div3)) / gamma_div
yangbase <- data.frame(place = c('固定沙地','半固定沙地','一年弃耕地'), shannon_winner = c(shannon1, shannon2, shannon3), simpson = c(simpson1, simpson2, simpson3), alpha_div = c(alpha_div1, alpha_div2, alpha_div3))
#label(yangbase$place) <- 'place'
#label(yangbase$shannon_winner) <- 'shannon-winner'
#label(yangbase$simpson) <- 'simpson'

#table1(~ place + shannon_winner + simpson ,dat=yangbase)

#空间格局

#多样性在空间上的变化
#x,y上
space_divs <- lapply(dyang, function(yid){
  ydat <- fildat(dat1, 1, yid)
  had <- levels(factor(ydat[,2]))
  alpha <- juan(matrix(1:100,10,10),function(x){
    x <- as.numeric(x)
    ans <- sapply(x, function(m){
      return(length(levels(factor(ydat[m,3]))))
    })
    return(mean(ans))
  })
  gamma <- juan(matrix(1:100,10,10),function(x){
    x <- as.numeric(x)
    return(length(levels(factor(ydat[x,3]))))
  })
  beta <- 1 - alpha[,2] / gamma[,2]
  beta <- ifelse(is.na(beta), 0, beta)
  return(cbind(alpha, gamma[,2], beta))
})

total_space_div <- space_divs[[1]]
for(each in 2:yangnum){
  total_space_div <- total_space_div + space_divs[[each]]
}
total_space_div <- total_space_div / yangnum
space_div <- sapply(levels(factor(total_space_div[,1])), function(space){
  tdat <- fildat(total_space_div, 1, space)
  if(space != darea){
    return(as.numeric(c(space, sapply(2:4, function(u){return(mean(tdat[,u]))}))))
  } else {
    return(as.numeric(tdat))
  }
})

jpeg('./scale_abu_in_desert/figure4.jpg', width=200*8, height=200*24)
par(mfcol=c(3,1))
plot(space_div[1,],space_div[2,],'o',col=4,xlab='area(m^2)',ylab='alpha diversity')
plot(space_div[1,],space_div[3,],'o',col=5,xlab='area(m^2)',ylab='gamma diversity')
plot(space_div[1,],space_div[4,],'o',col=6,xlab='area(m^2)',ylab='beta diversity')
dev.off()

#种-相对多度
jpeg('./scale_abu_in_desert/figure5.jpg', width=200*12, height=200*4)
par(mfcol=c(1,3))
m <- seq(0.1,1,0.1)
n <- sapply(seq(0.1,1,0.1), function(x){
  return(length(which(pis1<=x & pis1 > x - 0.1)))
})
plot(m, n, main = 'spsnum-abu-place1', xlab = 'abu rank', ylab = 'sps num', col = 4)
w <- coef(nls(n~a*m^b,data=as.data.frame(cbind(m,n))))
lines(m, sapply(m, function(x){
  return(w[1]*x^w[2])
}), lty = 2)

n <- sapply(seq(0.1,1,0.1), function(x){
  return(length(which(pis2<=x & pis2 > x - 0.1)))
})
plot(m, n, main = 'spsnum-abu-place2', xlab = 'abu rank', ylab = 'sps num', col = 5)
w <- coef(nls(n~a*m^b,data=as.data.frame(cbind(m,n))))
lines(seq(0.05,1,0.1), sapply(m, function(x){
  return(w[1]*x^w[2])
}), lty = 2)

n <- sapply(seq(0.1,1,0.1), function(x){
  return(length(which(pis3<=x & pis3 > x - 0.1)))
})
plot(m, n, main = 'spsnum-abu-place3', xlab = 'abu rank', ylab = 'sps num', col = 6)
w <- coef(nls(n~a*m^b,data=as.data.frame(cbind(m,n))))
lines(seq(0.05,1,0.1), sapply(m, function(x){
  return(w[1]*x^w[2])
}), lty = 2)

dev.off()


##植物大小和多度关系
#一个群落内，植物个体数目在植物个体大小上的分布
v_step <- 0.1

jpeg('./scale_abu_in_desert/figure6.jpg', width=200*8, height=200*4)
par(mfcol=c(1,3))
#b c spv1 <- scale(as.numeric(apply(dat1, 1, calV)))
spv1 <- log(as.numeric(apply(dat1, 1, calV)))
plot(seq(min(spv1),max(spv1),v_step), sapply(seq(min(spv1),max(spv1),v_step), function(x){
  return(length(which(spv1<=x & spv1 > x - qun_step)))
}), main = 'v-ind num-place1', xlab = 'log(V)', ylab = 'ind num', col = 5)

#spv2 <- scale(as.numeric(apply(dat2, 1, calV)))
spv2 <- log(as.numeric(apply(dat2, 1, calV)))
plot(seq(min(spv2),max(spv2),v_step), sapply(seq(min(spv2),max(spv2),v_step), function(x){
  return(length(which(spv2<=x & spv2 > x - qun_step)))
}), main = 'v-ind num-place2', xlab = 'log(V)', ylab = 'ind num', col = 6)

#spv3 <- scale(as.numeric(apply(dat3, 1, calV)))
spv3 <- log(as.numeric(apply(dat3, 1, calV)))
plot(seq(min(spv3),max(spv3),v_step), sapply(seq(min(spv3),max(spv3),v_step), function(x){
  return(length(which(spv3<=x & spv3 > x - qun_step)))
}), main = 'v-ind num-place3', xlab = 'log(V)', ylab = 'ind num', col = 7)

dev.off()

#逻辑斯蒂拟合
jpeg('./scale_abu_in_desert/figure7.jpg', width=200*8, height=200*4)
par(mfcol=c(1,3))

logis1 <- sapply(1:(length(levels(factor(sps1)))), function(n){
  spn <- sps1[n]
  spdat <- fildat(dat1, 3, spn)
  spv <- log(as.numeric(apply(spdat, 1, calV)))
  if(length(spv)<5){
    return(0)
  }
  m <- seq(min(spv),max(spv),qun_step)
  n <- sapply(seq(min(spv),max(spv),qun_step), function(x){
    return(length(which(spv<=x)))
  })
  ans <- sapply(2:length(n), function(x){return((n[x] - n[x-1]) / qun_step)})
  ans <- rev(ans)
  ans <- mean(ans[1:10])
  names(ans) <- spn
  return(ans)
})
plot(logis1[which(logis1 != 0)], pis1[which(logis1 != 0)], xlab = 'D', ylab = 'sp abu', main = paste('p-value',cor.test(logis1, pis1)$p.value), col = 4, pch = 15)
w <- my_lm(logis1[which(logis1 != 0)], pis1[which(logis1 != 0)])
lines(seq(min(logis1[which(logis1 != 0)]),max(logis1[which(logis1 != 0)]), 1), sapply(seq(min(logis1[which(logis1 != 0)]),max(logis1[which(logis1 != 0)]), 1), function(x){
  return(w[1] + w[2]*x)
}), lty = 2)

logis2 <- sapply(1:(length(levels(factor(sps2)))), function(n){
  spn <- sps2[n]
  spdat <- fildat(dat2, 3, spn)
  spv <- log(as.numeric(apply(spdat, 1, calV)))
  if(length(spv)<5){
    return(0)
  }
  m <- seq(min(spv),max(spv),qun_step)
  n <- sapply(seq(min(spv),max(spv),qun_step), function(x){
    return(length(which(spv<=x)))
  })
  ans <- sapply(2:length(n), function(x){return((n[x] - n[x-1]) / qun_step)})
  ans <- rev(ans)
  ans <- mean(ans[1:10])
  names(ans) <- spn
  return(ans)
})
plot(logis2[which(logis2 != 0)], pis2[which(logis2 != 0)], xlab = 'D', ylab = 'sp abu', main = paste('p-value',cor.test(logis2, pis2)$p.value,'*'), col = 3, pch = 15)
w <- my_lm(logis2[which(logis2 != 0)], pis2[which(logis2 != 0)])
lines(seq(min(logis2[which(logis2 != 0)]),max(logis2[which(logis2 != 0)]), 1), sapply(seq(min(logis2[which(logis2 != 0)]),max(logis2[which(logis2 != 0)]), 1), function(x){
  return(w[1] + w[2]*x)
}), lty = 2)

logis3 <- sapply(1:(length(levels(factor(sps3)))), function(n){
  spn <- sps3[n]
  spdat <- fildat(dat3, 3, spn)
  spv <- log(as.numeric(apply(spdat, 1, calV)))
  if(length(spv)<5){
    return(0)
  }
  m <- seq(min(spv),max(spv),qun_step)
  n <- sapply(seq(min(spv),max(spv),qun_step), function(x){
    return(length(which(spv<=x)))
  })
  ans <- sapply(2:length(n), function(x){return((n[x] - n[x-1]) / qun_step)})
  ans <- rev(ans)
  ans <- mean(ans[1:10])
  names(ans) <- spn
  return(ans)
})
tx <- logis3[which(logis3 != 0)] 
ty <- pis3[which(logis3 != 0)]
plot(tx, ty, xlab = 'D', ylab = 'sp abu', main = paste('p-value',cor.test(logis3, pis3)$p.value,'*'), col = 5 , pch = 15)
w <- my_lm(tx, ty)
lines(seq(min(tx),max(tx), 1), sapply(seq(min(tx),max(tx), 1), function(x){
  return(w[1] + w[2]*x)
}), lty = 2)
dev.off()

#植物个体大小-多度
mspv1 <- sapply(sps1, function(sp){
  return(mean(spv1[which(dat1[,3]==sp)]))
})
w <- my_lm(mspv1, pis1)
plot(mspv1, pis1, xlab = 'mean log(v)', ylab = 'sp abu', main = paste('p-value',cor.test(mspv1, pis1)$p.value), col = 5)
lines(seq(min(mspv1),max(mspv1), 1), sapply(seq(min(mspv1),max(mspv1), 1), function(x){
  return(w[1] + w[2]*x)
}), lty = 2)

mspv2 <- sapply(sps2, function(sp){
  return(mean(spv2[which(dat2[,3]==sp)]))
})
w <- my_lm(mspv2, pis2)
plot(mspv2, pis2, xlab = 'mean log(v)', ylab = 'sp abu', main = paste('p-value',cor.test(mspv2, pis2)$p.value), col = 5)
lines(seq(min(mspv2),max(mspv2), 1), sapply(seq(min(mspv2),max(mspv2), 1), function(x){
  return(w[1] + w[2]*x)
}), lty = 2)

mspv3 <- sapply(sps3, function(sp){
  return(mean(spv3[which(dat3[,3]==sp)]))
})
w <- my_lm(mspv3, pis3)
plot(mspv3, pis3, xlab = 'mean log(v)', ylab = 'sp abu', main = paste('p-value',cor.test(mspv3, pis3)$p.value), col = 5)
lines(seq(min(mspv3),max(mspv3), 1), sapply(seq(min(mspv3),max(mspv3), 1), function(x){
  return(w[1] + w[2]*x)
}), lty = 2)

#个体生物量和多度关系
biosp <- levels(factor(biodat[,1]))
biodat$cw_max <- ifelse(biodat$cw_max==0,1,biodat$cw_max)
biodat$cw_min <- ifelse(biodat$cw_min==0,1,biodat$cw_min)
biodat$height <- ifelse(biodat$height==0,1,biodat$height)
biodat$v <- log((biodat$height * biodat$cw_max * biodat$cw_min) / 4 * pi)
plot(0,0,'n',xlim=c(min(biodat$v),max(biodat$v)),ylim=c(min(biodat$biomass),max(biodat$biomass)), xlab = 'log(V)', ylab = 'biomass(g)')
sapply(1:length(biosp), function(sp){
  spn <- biosp[sp]
  pdat <- fildat(biodat, 1, spn)
  lines(pdat$v, pdat$biomass, 'p',col = sp)
})
legend('topleft',legend = biosp, col=1:length(biosp), pch=16)
w <- my_lm(biodat$v, biodat$biomass, 3)
lines(seq(0,12,0.5),sapply(seq(0,12,0.5),function(x){return(w[1]+w[2]*x+w[3]*x^2+w[4]*x^3)}))




