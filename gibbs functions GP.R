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
  Matrix(dmat)
  
#   teste=apply(dmat,2,sum)
#   plot(teste,table(dat$loc.id))
#   unique(teste-table(dat$loc.id))
}
#--------------------------
create.K=function(theta,covmat,nclust,nparam){
  tmp=matrix(0,nclust,nclust)
  for (i in 1:nparam){
    tmp=tmp+theta[i]*covmat[[i]]
  }
  exp(-tmp)
}
#--------------------------
update.alpha=function(param,dtd,xmat,dmat){
  prec=dtd+(1/param$sig2)*param$invK
  var1=as.matrix(solve(prec))
  err=param$z-xmat%*%param$betas
  pmedia=t(dmat)%*%err
  if (!isSymmetric(var1)) var1=(var1+t(var1))/2
  t(rmvnorm(1,var1%*%pmedia,var1))
}
#--------------------------
update.betas=function(param,xtx,xmat,dat){
  prec=xtx+diag(1,ncol(xmat))
  var1=solve(prec)
  err=param$z-param$alpha[dat$loc.id]
  pmedia=t(xmat)%*%err
  t(rmvnorm(1,var1%*%pmedia,var1))
}
#--------------------------
update.sig2=function(param,nclust,a.sig2,b.sig2){
  a=(nclust+2*a.sig2)/2
  err=param$alpha
  b=b.sig2+(t(err)%*%param$invK%*%err/2)  
  1/rgamma(1,a,b)
}
#--------------------------
update.z=function(param,dat,xmat){
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
get.inverse=function(D,sig2,K,dtd,nclust){
  med=solve(diag(1,nclust)+sig2*K%*%dtd)
  zzz=-sig2*D%*%med%*%K%*%t(D)
  diag(zzz)=1+diag(zzz)
  zzz
}
#--------------------------
update.theta=function(param,jump,xmat,nparam,ntheta.vals,nclust,dmat,dtd,theta.vals,covmat){
  theta.orig=theta.new=theta.old=param$theta
  err=param$z-xmat%*%param$betas
  
  for (i in 1:nparam){
    ind=which(theta.old[i]==theta.vals)
    seq1=(-jump[i]):jump[i]
    ind1=ind+sample(seq1,size=1)

    #make reflection
    if (ind1>ntheta.vals) {e=ind1-ntheta.vals; ind1=ntheta.vals-e}
    if (ind1<1)           {e=1-ind1;           ind1=1+e}

    theta.new[i]=theta.vals[ind1]
  }
  
  #calculate stuff
  K.old=param$K
  invK.old=param$invK
  K.new=create.K(theta.new,covmat,nclust,nparam)
  invK.new=solve(K.new)

  #---------------------
  tmp=(1/param$sig2)*invK.old+dtd
  p1.old=determinant(tmp,logarithm = T)$modulus[[1]]+
         determinant(param$sig2*K.old,logarithm = T)$modulus[[1]]
  inv.old=get.inverse(dmat,param$sig2,K.old,dtd,nclust)
  
  tmp=(1/param$sig2)*invK.new+dtd
  p1.new=determinant(tmp,logarithm = T)$modulus[[1]]+
         determinant(param$sig2*K.new,logarithm = T)$modulus[[1]]
  inv.new=get.inverse(dmat,param$sig2,K.new,dtd,nclust)
  
  # zzz=Matrix(diag(1,nobs))+param$sig2*dmat%*%K.new%*%t(dmat)
  # zzz1=solve(zzz)
  # hist(data.matrix(inv.new)-data.matrix(zzz1))
  # determinant(zzz,logarithm = T)$modulus[[1]]
  
  p2.old=t(err)%*%inv.old%*%err
  p2.new=t(err)%*%inv.new%*%err
  #---------------------
  pold=-(1/2)*(p1.old+p2.old)
  pnew=-(1/2)*(p1.new+p2.new)

  k=acceptMH(pold,pnew,theta.old,theta.new,T)
  theta.old=k$x

  if (k$accept==0) {K=K.old; invK=invK.old}
  if (k$accept==1) {K=K.new; invK=invK.new}

  list(theta=theta.old,K=K,invK=invK,accept=rep(k$accept,nparam))
}
