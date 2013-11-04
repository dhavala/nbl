library(mvtnorm)
library(moments)
library(VGAM)

# nbl random numbers
rnbl <- function(n,r,theta)
{
  
    lambda <- rlind(n,theta)
    p <- exp(-lambda)
    x <- rnbinom(n,r,prob=p)
}

# mgf of lindley distribution
mgflind <- function(x,theta)
{
  mgf <- (theta-x+1)/(theta-x)^2
  mgf <- mgf * (theta^2)/(theta+1)
}

# pmf at given vector
dnbl <- function(x,r,theta)
{
  
  sq <- seq(0,max(x))
  nx <- length(x)
  pr <- rep(0,nx)
  # for 0 to max(X), get the mgf
  mgf.seq <- mgflind(-(sq+r),theta)
  sign.seq <- (-1)^sq  
  for (ii in seq(1,nx))
  {
    k.local <- x[ii]
    ind <- c(0:k.local)
    nck.seq <- choose(k.local,ind)
    ip.seq <- sum(nck.seq*mgf.seq[ind+1]*sign.seq[ind+1])
    pr[ii] <- choose(r+k.local-1,k.local)*ip.seq
    
  }
  return(pr)
}

# pmf at a point
dnbl2 <- function(x,r,theta)
{
  
  jj <- seq(0,x)
  jjr <- r+jj
  yck <- choose(x,jj)
  signj <- (-1)^jj
  mgf <- (theta+jjr+1)*(theta/((theta+jjr)))^2/(theta+1)
  pr <- choose(r+x-1,x)*sum(yck*signj*mgf)
  return(pr)
}

# cdf
pnbl <- function(q,r,theta)
{
  
  jj <- seq(0,max(q))
  pr <- dnbl(jj,r,theta)
  cr <- cumsum(pr)
  return(cr[q+1])
}

# inverse-cdf
qnbl <- function(p,r,theta)
{
  # search within 15
  jj <- 15
  n <- length(p)
  q <- rep(NA, n)
  pvec <- pnbl(seq(0,jj),r,theta)
  if(pvec[jj+1]<max(p))
  {
    print('inverse cdf is inefficient')
    
  } else {
    
    for(ii in seq(1,n))
    {
      q[ii] <- which(pvec>=p[ii])[1]
    }  
  }
  return(q)
}

# example
theta <- 100
r <- 20
n = 1000

x <- rnbl(n,r,theta)
xu <- sort(uniq(x))
h <- hist(x,breaks=xu,plot=F)
h$counts <- h$counts/sum(h$counts)
plot(h)

yu <- dnbl(xu,r,theta)
lines(xu,yu,lwd=2,col='red')
mn <- mean(xu*yu)
lines(xu,dpois(xu,lambda=mn),lwd=2,col='blue')


dnbl(c(0,1,2),r,theta)
pnbl(c(0,1,2),r,theta)
qnbl(c(0.1,0.5,0.88,0.99),r,theta)