rm(list=ls(all=TRUE))
set.seed(18)
library('mvtnorm')
library("Matrix")

# setwd('U:\\ghana\\fever and malaria\\single DHS survey')
# dat=read.csv('actual data.csv',as.is=T)

setwd('U:\\independent studies\\gaussian process\\github_GP')
source('gibbs functions GP.R')
source('gibbs wrapper.R')
dat=read.csv('fake data.csv',as.is=T)

nomes.cov.vil=nomes=c('nl','elevation','mean_evi','dist_water','dist_urb','LATNUM1','LONGNUM1')
nomes.cov.ind='agekr.yr'
ngibbs=100
nburnin=10
res=GP.gibbs(dat,nomes.cov.vil,nomes.cov.ind,ngibbs,nburnin,include.pred=F)
  
vec.alpha=res$alpha
vec.theta=res$theta
vec.sig2=res$sig2
vec.betas=res$betas