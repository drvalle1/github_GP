compare1=function(true,estim){
  rango=range(c(true,estim))
  plot(true,estim,xlim=rango,ylim=rango)
  lines(rango,rango)
}

seq1=1000:1400
compare1(alpha.true,colMeans(vec.alpha[seq1,]))
compare1(theta.true,colMeans(vec.theta[seq1,]))
compare1(z.true,param$z)

#look at sig2: 
plot(vec.sig2[seq1],type='l')
abline(h=sigma2.true,col='red')

plot(density(vec.sig2[seq1]),type='l')
abline(v=sigma2.true)

#look at beta:
plot(vec.betas[seq1],type='l')
abline(h=betas.true,col='red')

#look more closely to theta:
# seq1=1:900
par(mfrow=c(3,3),mar=c(4,4,1,1))
for (i in 1:ncol(vec.theta)) {
  x=vec.theta[seq1,i]
  plot(x,type='l') #,ylim=quantile(x,c(0.01,0.99))
}
cor(vec.theta[seq1,])

theta.med=apply(vec.theta[seq1,],2,mean,na.rm=T)
par(mfrow=c(1,1),mar=rep(4,4))
rango=range(c(theta.true,theta.med))
plot(theta.true,theta.med,xlim=rango,ylim=rango)
lines(rango,rango)

# Sigma.new=create.K(theta.med,covmat,nclust,nparam)
# Sigma.true=create.K(theta.true,covmat,nclust,nparam)
# plot(Sigma.true,Sigma.new)
