# library(mvtnorm)
library(RobustGaSP)
library(RobustCalibration)
library(ggplot2)

num_obs=200

input=as.vector(seq(0,1,1/(num_obs-1)))

sigma_2=1 ##the result won't change no matter you change sigma_2


lambda_z=1000

R0=abs(outer(input,input,'-'))
eigenvalue_GaSP=matrix(NA,num_obs,3)
eigenvalue_SGaSP=matrix(NA,num_obs,3)
ecdf_GaSP_list=as.list(1:3)
ecdf_SGaSP_list=as.list(1:3)

for(i in 0:2 ){
  beta=.5*(10^i)
  
  R=matern_5_2_funct(R0,beta_i=beta)
  eigen_R=eigen(R)

  eigenvalue_GaSP[,i+1]=eigen_R$values/num_obs
  eigenvalue_SGaSP[,i+1]=eigen_R$values/(1+lambda_z/(num_obs)*eigen_R$values)/(num_obs)
  ###
  R_z=R-R%*%solve(R+num_obs/lambda_z*diag(num_obs))%*%R
  eigen_R_z=eigen(R_z)
  
  L=t(chol(R))
  L_z=t(chol(R_z))

  
  ratio=sum(eigen_R$values)/sum(eigen_R_z$values)
  
  M=500 ##number of simlation
  norm_sample=matrix(rnorm(M*num_obs),num_obs,M)
  output_GaSP=L%*%norm_sample
  Z_GaSP=colMeans(output_GaSP^2)
  
  output_SGaSP=L_z%*%norm_sample
  Z_SGaSP=colMeans(output_SGaSP^2*ratio)
  
  ecdf_GaSP_list[[i+1]]=ecdf(Z_GaSP)
  ecdf_SGaSP_list[[i+1]]=ecdf(Z_SGaSP)
  
}

##note these two are the same
eigen_R$values/(1+lambda_z/(num_obs)*eigen_R$values)/(num_obs)-eigen_R_z$values/(num_obs)


 pdf('eigenvalues_comparison.pdf',height=4,width=5)
  plot(log(eigenvalue_GaSP[,1]),col='red',pch=0,ylab=expression(log(eigenvalue[i])),xlab='i',mgp=c(2.5,1,0))
  
  lines(log(eigenvalue_GaSP[,2]),col='red',pch=1,type='p')
    lines(log(eigenvalue_SGaSP[,1]),col='blue',pch=15,type='p')
  lines(log(eigenvalue_SGaSP[,2]),col='blue',pch=19,type='p')
  legend('topright',legend=c(expression(GaSP~~gamma==2),
                             expression(GaSP~~gamma==.2),expression(S-GaSP~~gamma==2),
                             expression(S-GaSP~~gamma==.2)),pch=c(0,1,15,19),col=c('red','red','blue','blue'))
  dev.off()

  pdf('Z_CDF_comparison.pdf',height=4,width=5)
  plot(ecdf_GaSP_list[[1]],col='red',pch=0,ylab='CDF',xlab='Z',mgp=c(2.5,1,0),main='')
  lines(ecdf_GaSP_list[[2]],col='red',pch=1)
  lines(ecdf_SGaSP_list[[1]],col='blue',pch=15)
  lines(ecdf_SGaSP_list[[2]],col='blue',pch=19)
  
  legend('bottomright',legend=c(expression(GaSP~~gamma==2),
                             expression(GaSP~~gamma==.2),expression(variance~adjusted~S-GaSP~~gamma==2),
                             expression(variance~adjusted~S-GaSP~~gamma==.2)),pch=c(0,1,15,19),col=c('red','red','blue','blue'))
  dev.off()
  

