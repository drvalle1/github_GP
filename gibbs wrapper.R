GP.gibbs=function(dat,nomes.cov.vil,nomes.cov.ind,ngibbs,nburnin,include.pred){
  nobs=nrow(dat)
  
  #for prediction
  ind.pred=which(is.na(dat$microsc1))
  npred=length(ind.pred)
  
  #get village data
  nomes=c('loc.id',nomes.cov.vil)
  clust=unique(round(dat[,sort(nomes)],4))
  clust=clust[order(clust$loc.id),]
  nclust=nrow(clust); nclust
  
  #create covariate "distance" matrices
  covmat=list()
  for (i in 1:length(nomes.cov.vil)){
    vec=clust[,nomes.cov.vil[i]]
    covmat[[i]]=outer(vec,vec,'-')^2
  }
  names(covmat)=nomes.cov.vil
  nparam=length(covmat)
  
  #get dmat
  dmat=create.dmat(dat)
  dtd=t(dmat)%*%dmat
  
  #get xmat and betas
  xmat=data.matrix(dat[,nomes.cov.ind])
  betas=rep(0,ncol(xmat))
  xtx=t(xmat)%*%xmat
  
  #discretize theta to improve mixing
  cors=seq(from=0.0001,to=0.9,length.out=11)
  res=matrix(0,nclust,nclust)
  for (i in 1:length(covmat)) res=res+covmat[[i]]
  tmp1=c(-log(cors)/quantile(res,0.5))
  tmp2=c(-log(cors)/quantile(res,0.1))
  theta.vals=sort(unique(c(tmp1,tmp2)))
  ntheta.vals=length(theta.vals)
  theta=rep(theta.vals[11],nparam)
  K=create.K(theta,covmat,nclust,nparam)
  invK=solve(K)
  
  #get initial values
  alpha=rep(0,nclust)
  sig2=2
  a.sig2=b.sig2=0.1
  
  #get z
  z=rep(NA,nrow(dat))
  cond=!is.na(dat$microsc1) & dat$microsc1;    z[cond]=1
  cond=!is.na(dat$microsc1) & dat$microsc1==0; z[cond]=-1
  cond=is.na(dat$microsc1);                    z[cond]=-0.5;
  
  #useful things for gibbs sampler
  vec.alpha=matrix(NA,ngibbs,nclust)
  vec.theta=matrix(NA,ngibbs,nparam)
  vec.sig2=matrix(NA,ngibbs,1) #sig2
  vec.betas=matrix(NA,ngibbs,length(betas))
  vec.pred=matrix(NA,ngibbs,npred)
  
  param=list(alpha=matrix(alpha,nclust,1),theta=theta,betas=betas,z=z,sig2=sig2,
             invK=invK,K=K)
  
  jump=rep(1,nparam)
  accept.theta=rep(0,nparam)
  accept.output=100
  
  options(warn=2)
  for (i in 1:ngibbs){
    print(i)
    # if (i==1){
    #   param$theta=theta.true
    #   param$Sigma=create.Sigma(theta.true,covmat,nclust)
    #   param$invSigma=solve(param$Sigma)
    # }
    tmp=update.theta(param,jump,xmat,nparam,ntheta.vals,nclust,dmat,dtd,theta.vals,covmat)
    param$theta=tmp$theta
    param$invK=tmp$invK
    param$K=tmp$K
    accept.theta=accept.theta+tmp$accept
    
    param$alpha=update.alpha(param,dtd,xmat,dmat)#t(alpha.true)#
    
    param$sig2=update.sig2(param,nclust,a.sig2,b.sig2) #sigma2.true#
    param$betas=update.betas(param,xtx,xmat,dat) #betas.true#
    
    #sample states
    param$z=update.z(param,dat,xmat)#z.true#
    
    if (i%%accept.output==0 & i<nburnin){
      rate=accept.theta/accept.output
      print(rate)
      print(jump)
      #   cond=rate>0.5; jump[cond]=jump[cond]+3
      #   cond=rate<0.1 & (jump-3)>0
      #   jump[cond]=jump[cond]-3
      accept.theta=rep(0,nparam)
      #   jump[]=1
    }
    
    #store results
    vec.alpha[i,]=param$alpha
    vec.theta[i,]=param$theta
    vec.sig2[i]=param$sig2 
    vec.betas[i,]=param$betas
    
    #for prediction
    if (include.pred){
      media=param$alpha[dat$loc.id]+xmat%*%param$betas
      vec.pred[i,]=pnorm(media[ind.pred])
    }
    
  }
  list(alpha=vec.alpha,theta=vec.theta,sig2=vec.sig2,betas=vec.betas,pred=vec.pred)  
}
