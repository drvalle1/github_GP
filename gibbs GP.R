# rm(list=ls(all=TRUE))
set.seed(18)
library('mvtnorm')
library("Matrix")

# setwd('U:\\ghana\\fever and malaria\\single DHS survey')
# dat=read.csv('actual data.csv',as.is=T)

setwd('U:\\independent studies\\gaussian process\\github_GP')
source('gibbs functions GP.R')
dat=read.csv('fake data.csv',as.is=T)
nobs=nrow(dat)

#get village data
tipos=''
nomes=c('LATNUM1','LONGNUM1','loc.id',
        paste('nl',tipos,sep=''),
        paste('elevation',tipos,sep=''),
        paste('mean_evi',tipos,sep=''),
        paste('dist_water',tipos,sep=''),
        paste('dist_urb',tipos,sep=''))
clust=unique(round(dat[,sort(nomes)],4))
clust=clust[order(clust$loc.id),]
nclust=nrow(clust); nclust

#create covariate "distance" matrices
nomes.cov=nomes=c('nl','elevation','mean_evi','dist_water','dist_urb','LATNUM1','LONGNUM1')
covmat=list()

for (i in 1:length(nomes.cov)){
  vec=clust[,nomes.cov[i]]
  covmat[[i]]=outer(vec,vec,'-')^2
}
names(covmat)=nomes.cov
nparam=length(covmat)

#get dmat
dmat=create.dmat(dat)
dtd=t(dmat)%*%dmat

#get xmat and betas
nomes.x=c('agekr.yr')
xmat=data.matrix(dat[,nomes.x])
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
ngibbs=10000
vec.alpha=matrix(NA,ngibbs,nclust)
vec.theta=matrix(NA,ngibbs,nparam)
vec.sig2=matrix(NA,ngibbs,1) #sig2
vec.betas=matrix(NA,ngibbs,length(betas))

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
  tmp=update.theta(param,jump)
  param$theta=tmp$theta
  param$invK=tmp$invK
  param$K=tmp$K
  accept.theta=accept.theta+tmp$accept

  param$alpha=update.alpha(param)#t(alpha.true)#
  
  param$sig2=update.sig2(param) #sigma2.true#
  param$betas=update.betas(param) #betas.true#
  
  #sample states
  param$z=update.z(param)#z.true#

  if (i%%accept.output==0 & i<1000){
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
}
