compare1=function(true,estim){
  rango=range(c(true,estim))
  plot(true,estim,xlim=rango,ylim=rango)
  lines(rango,rango)
}

seq1=1000:9900
compare1(alpha.true,colMeans(vec.alpha[seq1,]))
compare1(theta.true,param$theta)
compare1(z.true,param$z)

#look at sig2: 
plot(vec.sig2[seq1],type='l')
abline(h=sigma2.true,col='red')

plot(density(vec.sig2[seq1]),type='l')
abline(v=sigma2.true)

#look at beta:
plot(vec.betas[seq1],type='l')
abline(h=betas.true,col='red')