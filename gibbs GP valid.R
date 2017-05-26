rm(list=ls(all=TRUE))
set.seed(18)
library('mvtnorm')
library("Matrix")

setwd('U:\\independent studies\\gaussian process\\github_GP')
source('gibbs functions GP.R')
source('gibbs wrapper.R')

nomes.cov.vil=nomes=c('nl','elevation','mean_evi','dist_water','dist_urb','LATNUM1','LONGNUM1')
nomes.cov.ind='agekr.yr'
ngibbs=5000
nburnin=1000

for (valid in 1:10){
  setwd('U:\\independent studies\\gaussian process\\validation')
  nome=paste('valid',valid,'a.csv',sep='')
  dat=read.csv(nome,as.is=T)
  
  res=GP.gibbs(dat,nomes.cov.vil,nomes.cov.ind,ngibbs,nburnin,include.pred=T)

  ind.pred=which(is.na(dat$microsc1))
  seq1=nburnin:ngibbs
  fim=data.frame(ind=dat$ind[ind.pred],
                 prev.GP=colMeans(res$pred[seq1,]))
  
  setwd('U:\\independent studies\\gaussian process\\validation\\results')
  nome=paste('valid',valid,'_GP.csv',sep='')
  write.csv(fim,nome,row.names=F)
}
