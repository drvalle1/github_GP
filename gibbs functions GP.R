tnorm <- function(n,lo,hi,mu,sig){   #generates truncated normal variates based on cumulative normal distribution
  #normal truncated lo and hi
  
  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))
  
  q1 <- pnorm(lo,mu,sig) #cumulative distribution
  q2 <- pnorm(hi,mu,sig) #cumulative distribution
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  z[z == -Inf]  <- lo[z == -Inf]
  z[z == Inf]   <- hi[z == Inf]
  z
}
#----------------------------------------------------------------------------------------------
acceptMH <- function(p0,p1,x0,x1,BLOCK){   #accept for M, M-H
  # if BLOCK, then accept as a block,
  # otherwise, accept individually
  
  nz           <- length(x0)  #no. to accept
  if(BLOCK) nz <- 1
  
  a    <- exp(p1 - p0)       #acceptance PR
  z    <- runif(nz,0,1)
  keep <- which(z < a)
  
  if(BLOCK & length(keep) > 0) x0 <- x1
  if(!BLOCK)                   x0[keep] <- x1[keep]           
  accept <- length(keep)        
  
  list(x = x0, accept = accept)
}
#-------------------------------
rmvnorm1=function (n, sigma, pre0.9_9994 = FALSE) 
{
  #   retval <- chol(sigma, pivot = TRUE)
  #   o <- order(attr(retval, "pivot"))
  #   retval <- retval[, o]
  s. <- svd(sigma)
  if (!all(s.$d >= -sqrt(.Machine$double.eps) * abs(s.$d[1]))) {
    warning("sigma is numerically not positive definite")
  }
  R = t(s.$v %*% (t(s.$u) * sqrt(s.$d)))
  retval <- matrix(rnorm(n * ncol(sigma)), nrow = n, byrow = !pre0.9_9994) %*% R
  retval
}
#----------------------------
create.dmat=function(dat){
  nloc=length(unique(dat$loc.id))
  dmat=matrix(0,nrow(dat),nloc)
  for (i in 1:nloc){
    cond=dat$loc.id==i
    dmat[cond,i]=1
  }
  dmat
  
#   teste=apply(dmat,2,sum)
#   plot(teste,table(dat$loc.id))
#   unique(teste-table(dat$loc.id))
}
#--------------------------
create.Sigma=function(theta,covmat,nclust){
  tmp=matrix(0,nclust,nclust)
  for (i in 1:nparam){
    tmp=tmp+(covmat[[i]]/theta[i])
  }
  exp(-tmp)
}
#--------------------------
update.alpha=function(param){
  prec=dtd+(1/param$sig2)*param$invSigma
  var1=solve(prec)
  err=param$z-xmat%*%param$betas
  pmedia=t(dmat)%*%err
  t(rmvnorm(1,var1%*%pmedia,var1))
}
#--------------------------
update.betas=function(param){
  prec=xtx+diag(1,ncol(xmat))
  var1=solve(prec)
  err=param$z-param$alpha[dat$loc.id]
  pmedia=t(xmat)%*%err
  t(rmvnorm(1,var1%*%pmedia,var1))
}
#--------------------------
update.sig2=function(param){
  a=(nclust+2*a.sig2)/2
  err=param$alpha
  b=b.sig2+(t(err)%*%param$invSigma%*%err/2)  
  1/rgamma(1,a,b)
}
#--------------------------
update.z=function(param){
  media=param$alpha[dat$loc.id]+xmat%*%param$betas

  z=rep(NA,nrow(dat))
  cond=!is.na(dat$microsc1) & dat$microsc1==1
  z[cond]=tnorm(sum(cond),lo=0,hi=Inf ,mu=media[cond],sig=1)
  cond=!is.na(dat$microsc1) & dat$microsc1==0
  z[cond]=tnorm(sum(cond),lo=-Inf,hi=0,mu=media[cond],sig=1)
  cond= is.na(dat$microsc1)
  z[cond]=rnorm(sum(cond),mean=media[cond],sd=1)
  z
}
#--------------------------
update.theta=function(param){
  theta.new=theta.old=param$theta

  #I need to add probabilities associated with the bounds
  p=rep(NA,nparam)
  for (i in 1:nparam){
    ind=which(theta.old[i]==theta.vals)
    delta=sample(c(-1,0,1),size=1)  #move one step
    ind1=ind+delta
    if (ind1>ntheta.vals) ind1=ntheta.vals
    if (ind1<1)           ind1=1
    
    #special cases
    old.new=new.old=1
    if (ind==ntheta.vals & ind1==ntheta.vals-1) {old.new=1/2; new.old=1/3}
    if (ind==ntheta.vals-1 & ind1==ntheta.vals) {old.new=1/3; new.old=1/2} 
    if (ind==1 & ind1==2)             {old.new=1/2; new.old=1/3}
    if (ind1==2 & ind==1)             {old.new=1/3; new.old=1/2}
    p[i]=new.old/old.new
    
    theta.new[i]=theta.vals[ind1]
  }
  Sigma.new=create.Sigma(theta.new,covmat,nclust)
  invSigma.new=solve(Sigma.new)
  
  p1.old=determinant(param$Sigma,logarithm = T)$modulus[[1]]
  p1.new=determinant(Sigma.new  ,logarithm = T)$modulus[[1]]
  
  med.old=t(param$alpha)%*%param$invSigma%*%param$alpha/2
  med.new=t(param$alpha)%*%invSigma.new%*%param$alpha/2
  
  p2.old=(nclust+2*a.sig2)*log(med.old+b.sig2)
  p2.new=(nclust+2*a.sig2)*log(med.new+b.sig2)
  pold=-(1/2)*(p1.old+p2.old)
  pnew=-(1/2)*(p1.new+p2.new)
  
  k=acceptMH(pold,pnew+sum(log(p)),theta.old,theta.new,T)  
  if (k$accept==0) fim=list(theta=param$theta,Sigma=param$Sigma,invSigma=param$invSigma,accept=0)
  if (k$accept==1) fim=list(theta=theta.new  ,Sigma=Sigma.new  ,invSigma=invSigma.new  ,accept=1)
  fim
}
  