# rm(list=ls(all=TRUE))
set.seed(18)
library('mvtnorm')

# setwd('U:\\ghana\\fever and malaria\\single DHS survey')
# dat=read.csv('actual data.csv',as.is=T)

setwd('U:\\independent studies\\gaussian process')
source('gibbs functions GP.R')
dat=read.csv('fake data.csv',as.is=T)

#get village data
tipos=''
nomes=c('LATNUM','LONGNUM','loc.id',
        paste('nl',tipos,sep=''),
        paste('elevation',tipos,sep=''),
        paste('mean_evi',tipos,sep=''),
        paste('dist_water',tipos,sep=''),
        paste('dist_urb',tipos,sep=''))
clust=unique(round(dat[,sort(nomes)],4))
clust=clust[order(clust$loc.id),]
nclust=nrow(clust); nclust

#create covariate "distance" matrices
dist1=as.matrix(dist(clust[,c('LATNUM','LONGNUM')]))
nomes.cov=nomes=c('nl','elevation','mean_evi','dist_water','dist_urb')
covmat=list()

for (i in 1:length(nomes.cov)){
  vec=clust[,nomes.cov[i]]
  covmat[[i]]=outer(vec,vec,'-')^2
}
covmat[[i+1]]=dist1
names(covmat)=c(nomes.cov,'distance')
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
cors=seq(from=0.0001,to=0.999,length.out=100)
theta.vals=-2/log(cors)
ntheta.vals=length(theta.vals)
theta=sample(theta.vals,size=nparam)

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

Sigma=create.Sigma(theta,covmat,nclust)
invSigma=solve(Sigma)
param=list(alpha=matrix(alpha,nclust,1),theta=theta,betas=betas,z=z,sig2=sig2,
           invSigma=invSigma,Sigma=Sigma)
accept.theta=0
options(warn=2)
for (i in 1:ngibbs){
  print(i)
  param$alpha=update.alpha(param)#t(alpha.true)#
  
  # if (i==1){
  #   param$theta=theta.true
  #   param$Sigma=create.Sigma(theta.true,covmat,nclust)
  #   param$invSigma=solve(param$Sigma)
  # }
  tmp=update.theta(param)
  param$theta=tmp$theta
  param$invSigma=tmp$invSigma
  param$Sigma=tmp$Sigma
  accept.theta=accept.theta+tmp$accept
  
  param$sig2=update.sig2(param) #sigma2.true#
  param$betas=update.betas(param) #betas.true#
  
  #sample states
  param$z=update.z(param)#z.true#

  #store results
  vec.alpha[i,]=param$alpha
  vec.theta[i,]=param$theta
  vec.sig2[i]=param$sig2 
  vec.betas[i,]=param$betas
}
accept.theta

