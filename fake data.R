rm(list=ls(all=TRUE))
set.seed(5)
library('mvtnorm')

#get data
setwd('U:\\ghana\\fever and malaria\\single DHS survey')
dat2=read.csv('actual data.csv',as.is=T)

#create village level matrix
tipos=''
nomes=c('LATNUM','LONGNUM','loc.id',
        paste('nl',tipos,sep=''),
        paste('elevation',tipos,sep=''),
        paste('mean_evi',tipos,sep=''),
        paste('dist_water',tipos,sep=''),
        paste('dist_urb',tipos,sep=''))
clust=unique(round(dat2[,sort(nomes)],4))
clust=clust[order(clust$loc.id),]; nrow(clust)

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

#create correlation matrices
cors=seq(from=0.0001,to=0.999,length.out=11)
theta.vals=-2/log(cors)
theta.true=theta=c(theta.vals[10],theta.vals[1],theta.vals[10],rep(theta.vals[1],3))
tmp=matrix(0,nrow(dist1),nrow(dist1))
for (i in 1:length(covmat)){
  tmp=tmp+(covmat[[i]]/theta[i])
}
Sigma=exp(-tmp)

#generate alphas
sigma2.true=sigma2=5
alpha.true=alpha=rmvnorm(1,rep(0,nrow(Sigma)),sigma2*Sigma)

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

setwd('U:\\independent studies\\gaussian process')
write.csv(dat3,'fake data.csv',row.names=F)
