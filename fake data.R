rm(list=ls(all=TRUE))
set.seed(3)
library('mvtnorm')

setwd('U:\\independent studies\\gaussian process\\github_GP')
source('gibbs functions GP.R')

#get data
setwd('U:\\ghana\\fever and malaria\\single DHS survey')
dat2=read.csv('actual data.csv',as.is=T)
dat2$LATNUM1=(dat2$LATNUM-mean(dat2$LATNUM))/sd(dat2$LATNUM)
dat2$LONGNUM1=(dat2$LONGNUM-mean(dat2$LONGNUM))/sd(dat2$LONGNUM)

#create village level matrix
tipos=''
nomes=c('LATNUM1','LONGNUM1','loc.id',
        paste('nl',tipos,sep=''),
        paste('elevation',tipos,sep=''),
        paste('mean_evi',tipos,sep=''),
        paste('dist_water',tipos,sep=''),
        paste('dist_urb',tipos,sep=''))
clust=unique(round(dat2[,sort(nomes)],4))
clust=clust[order(clust$loc.id),]; 
nclust=nrow(clust); nclust

#create covariate "distance" matrices
nomes.cov=nomes=c('nl','elevation','mean_evi','dist_water','dist_urb','LATNUM1','LONGNUM1')
covmat=list()

for (i in 1:length(nomes.cov)){
  vec=clust[,nomes.cov[i]]
  covmat[[i]]=outer(vec,vec,'-')^2
}
names(covmat)=nomes.cov

#create correlation matrices
cors=seq(from=0.0001,to=0.9,length.out=11)
res=matrix(0,nclust,nclust)
for (i in 1:length(covmat)) res=res+covmat[[i]]

theta.vals=c(-log(cors)/quantile(res,0.5))
tipos=c(rep(11,3),rep(1,length(covmat)-3))
ind=sample(tipos,size=length(covmat),replace=F)
theta.true=theta=theta.vals[ind]
Sigma=create.K(theta,covmat,nclust,nparam=length(covmat))
tmp=Sigma; diag(tmp)=NA
hist(tmp); mean(tmp>0.5,na.rm=T)

#generate alphas
sigma2.true=sigma2=3
alpha.true=alpha=rmvnorm(1,rep(0,nrow(Sigma)),sigma2*Sigma)
range(alpha.true)

#generate infection status
dat3=dat2
xmat=data.matrix(dat3[,'agekr.yr'])
betas.true=betas=-1#runif(ncol(xmat),min=-2,max=2)
media=alpha[dat3$loc.id]+xmat%*%betas
range(xmat%*%betas)
z.true=z=rnorm(nrow(dat3),mean=media,sd=1)
range(z)

#replace infected
mean(dat3$microsc1,na.rm=T)
dat3$microsc1=ifelse(z>0,1,0)
mean(dat3$microsc1,na.rm=T)

setwd('U:\\independent studies\\gaussian process\\github_GP')
write.csv(dat3,'fake data.csv',row.names=F)
