##rewrite this code by the  code from FastGaSP package and a new function for linear coefficient
#library(FastGaSP)
 library(Rcpp)
 library(RcppEigen)
# library(dlm)
# library(mvtnorm)
 library(RobustGaSP)
 sourceCpp(file='src/functions_June_30_2020.cpp') 

 
 y_R<-function(x){
   j_seq=seq(1,100,1)
   record_y_R=0
   for(i_j in 1:100){
     record_y_R=record_y_R+j_seq[i_j]^{-2*3}*sin(j_seq[i_j]*5)*cos(pi*(j_seq[i_j]-0.5)*x*5)

   }
   record_y_R
 }
 

pred_mean_FFBS_X_GP<-function(sigma_2,beta,eta, delta_x,index_obs,output){
  
  param=log(c(beta,eta))
  have_noise=T
  theta_hat=Get_linear_coeff_1_dim_GP(param, have_noise, delta_x, output,
                                    kernel_type,  X)

  return_list=as.list(1:2)
  
  return_list[[1]]=theta_hat
  return_list[[2]]=X_testing*theta_hat+Kalman_smoother(param, have_noise, index_obs, 
                                                 delta_x_all,  output-theta_hat, 1.0,
                                                 kernel_type,only_mean)[[1]]  
  return_list
}


pred_mean_FFBS_X_SGP<-function(sigma_2,beta,eta, lambda_z,delta_x,index_obs,output){
  
  param=log(c(beta,eta))
  have_noise=T
  theta_hat=Get_linear_coeff_1_dim_SGP(param, lambda_z,have_noise, delta_x, as.vector(output),
                                   kernel_type,  X)
  
  

  return_list=as.list(1:2)
  
  param2=log(c(beta,eta/shrinkage))
  
  
  return_list[[1]]=theta_hat
  

  return_list[[2]]=X_testing*theta_hat+1/shrinkage*Kalman_smoother(param2, have_noise, index_obs, 
                                                       delta_x_all,  output-theta_hat, 1.0,
                                                       kernel_type,only_mean)[[1]]  
  return_list
}


plot(y_R(seq(0,1,0.001)))

mean(y_R(seq(0,1,0.001)))


n_sim=50
n_experiment=50


n_record_theta_GP=matrix(0,n_sim,n_experiment)
n_record_theta_SGP=matrix(0,n_sim,n_experiment)
n_record_RMSE_GP=matrix(0,n_sim,n_experiment)
n_record_RMSE_SGP=matrix(0,n_sim,n_experiment)

delta_hat_L2=matrix(0,n_sim,n_experiment)
record_SGP_pred_sum=matrix(0,n_sim,n_experiment)

kernel_type='matern_5_2'
only_mean=T
#n_record_theta=matrix(0,n_sim,2)
#n_record_RMSE=matrix(0,n_sim,2)

num_obs_all=3*10^4


system.time(
  for(i_n in 1:n_sim){
    print(i_n)
    #set.seed(i_n)
    
    num_obs=floor(exp(5+i_n/10))
    #  num_obs*5
    #num_obs=100
    output_all=rep(NA,num_obs_all)
    
    index_obs=rep(0,num_obs_all)  ##zero to 1
    
    obs_index=floor(seq(1:num_obs)*num_obs_all/num_obs)
    
    input_all=seq(0,1,1/(num_obs_all-1))
    
    index_obs[obs_index]=1
    index_obs=as.integer(index_obs)
    
    input=input_all[obs_index]
    X=rep(1,num_obs)
    X_testing=rep(1,num_obs_all)
    ##########################here including the ones that need to be interpolate
    delta_x_all=input_all[2:num_obs_all]-input_all[1:(num_obs_all-1)];
    delta_x=input[2:num_obs]-input[1:(num_obs-1)];
    
    for(i_exp in 1:n_experiment){
      #print(i_exp)
      
      set.seed(i_exp)
      output=y_R(input)+rnorm(num_obs,mean=0,sd=0.05)
      #output_all[obs_index]=output
      testing_output=y_R(input_all)
      


      ##
      beta=1
      sigma_2=1 ##the result won't change no matter you change sigma_2
      #lambda=num_obs^{-1}*1*10^{-4}  ###this is not optimal...
      
      lambda=num_obs^{-2*3/(2*3+1)}*1*10^{-8}  ###this is not optimal...
      
      
      eta=num_obs*lambda
      
      lambda_z=sqrt(1/lambda)
      #lambda_z=(1/lambda)  ##this won't converge 
      
      shrinkage=1+lambda*lambda_z
      
      
      
      pred_mean_X_all_GP=pred_mean_FFBS_X_GP(sigma_2,beta,eta, delta_x,index_obs,output)
      pred_mean_X_all_SGP=pred_mean_FFBS_X_SGP(sigma_2,beta,eta,lambda_z, delta_x,index_obs,output)
      
      #plot(pred_mean_X_all_GP[[2]])
      #lines(testing_output,col='blue')
      
      #pred_mean_X_all_SGP=pred_mean_FFBS_X_SGP(sigma_2,beta,eta,lambda,lambda_z, delta_x,output_all,obs_index,output,theta)
      n_record_theta_GP[i_n,i_exp]=pred_mean_X_all_GP[[1]]  ##theta
      n_record_RMSE_GP[i_n,i_exp]=sqrt(mean((pred_mean_X_all_GP[[2]]-testing_output)^2))
      
      n_record_theta_SGP[i_n,i_exp]=pred_mean_X_all_SGP[[1]]  ##theta
      n_record_RMSE_SGP[i_n,i_exp]=sqrt(mean((pred_mean_X_all_SGP[[2]]-testing_output)^2))
      

    }
    print(c( mean(n_record_RMSE_GP[i_n,]),mean(n_record_RMSE_SGP[i_n,]) ))

  }
)



n_record_RMSE_GP_mean=matrix(0,n_sim,n_experiment )
n_record_RMSE_SGP_mean=matrix(0,n_sim,n_experiment )

for(i_sim in 1:n_sim){
  for(i_experiment in 1:n_experiment){
    n_record_RMSE_GP_mean[i_sim,i_experiment]=sqrt( mean((testing_output-n_record_theta_GP[i_sim,i_experiment])^2))
    n_record_RMSE_SGP_mean[i_sim,i_experiment]=sqrt( mean((testing_output-n_record_theta_SGP[i_sim,i_experiment])^2))
    
  }
}


library(ggplot2)
n_seq=exp(5+(1:n_sim)/10)


L_2_minimizer= mean(y_R(seq(0,1,10^{-5})))

index_set=1:50


record_RMSE_est_GP_SGP=data.frame(method = factor(rep(c("GaSP", "S-GaSP"), each=n_sim)), 
                                   RMSE = c(as.vector(rowMeans(n_record_RMSE_GP) ),as.vector(rowMeans(n_record_RMSE_SGP))))
record_RMSE_theory=data.frame(method = factor(rep(c("GaSP", "S-GaSP"), each=n_sim)), 
                              RMSE=c(as.vector((n_seq^{-3/(2*3+1)}/5) ),as.vector((n_seq^{-3/(2*3+1)}/5))))

pdf("pred_fixed_lambda_GP_SGP.pdf",heigh=3,width=4.5)
ggplot(record_RMSE_est_GP_SGP, aes(x=rep(log(n_seq),2), y=RMSE,shape=method,color=method))+
  labs(x ='log(n)',y=expression(AvgRMSE[f^M+delta]))+
  geom_point()+scale_colour_manual(name = "method", values = c("red", "blue"), 
                                                      labels = c('GaSP', 'S-GaSP')) +
geom_line(color='black',data = record_RMSE_theory, aes(x=rep(log(n_seq),2), y=RMSE),linetype="dashed")
dev.off()

record_theta_est_GP_SGP=data.frame(method = factor(rep(c("GaSP", "S-GaSP"), each=n_sim)), 
                                  RMSE = c(as.vector(log(sqrt(rowMeans((n_record_theta_GP-L_2_minimizer)^2) ))),as.vector(log(sqrt(rowMeans((n_record_theta_SGP-L_2_minimizer)^2)) ))))
record_RMSE_theta_theory=data.frame(method = factor(rep(c("GaSP", "S-GaSP"), each=n_sim)), 
                              RMSE=c(as.vector(log((n_seq^{-3/(2*3+1)}/40)) ),as.vector( log((n_seq^{-3/(2*3+1)}/40)))))

pdf("theta_fixed_lambda_GP_SGP.pdf",heigh=3,width=4.5)

ggplot(record_theta_est_GP_SGP, aes(x=rep(log(n_seq),2), y=RMSE,shape=method,color=method))+
  labs(x ='log(n)',y=expression(log(RMSE[theta])))+
  geom_point()+scale_colour_manual(name = "method", values = c("red", "blue"), 
                                   labels = c('GaSP', 'S-GaSP')) +
geom_line(color='black',data = record_RMSE_theta_theory, aes(x=rep(log(n_seq),2), y=RMSE),linetype="dashed")
dev.off()











pdf("pred_fixed_lambda_GP_SGP.pdf",heigh=5,width=7)
plot(log(n_seq)[index_set],(rowMeans(n_record_RMSE_GP)[index_set]),type='p',col='red',ylim=c(0,0.03),xlab=expression(log(n)),ylab=expression(AvgRMSE[f^M+delta]),pch=2,cex=0.8)
lines(log(n_seq)[index_set],(rowMeans(n_record_RMSE_SGP)[index_set]),type='p',col='blue',pch=1,cex=0.8)
##m/(2m+1) m=2.5
#lines(log(n_seq),n_seq^{-2.5/(2*2.5+1)}/10.5,type='l')
lines(log(n_seq)[index_set],((n_seq^{-3/(2*3+1)}/4.9)[index_set]),type='l',cex=0.8)
legend("bottomleft", 
       legend = c("GaSP", "S-GaSP"), 
       col = c('red','blue'),
       pch=c(2,1),cex=0.8)
dev.off()


pdf("theta_fixed_lambda_GP_SGP.pdf",heigh=5,width=7)
plot(log(n_seq),log(sqrt(rowMeans(n_record_theta_GP-L_2_minimizer)^2))[index_set],type='p',col='red',ylim=c(-8,5),
     ,xlab=expression(log(n)),ylab=expression(log(RMSE[theta])),pch=2,cex=0.8)
lines(log(n_seq)[index_set],log(sqrt(rowMeans((n_record_theta_SGP-L_2_minimizer)^2)))[index_set],type='p',col='blue',
      pch=1,cex=0.8)
# plot(log(n_seq),log(rowMeans(abs(n_record_theta_GP-L_2_minimizer))),type='p',col='red',ylim=c(-8,1),
#      ,xlab=expression(log(n)),ylab=expression(log(RMSE[theta])),pch=1,cex=0.8)
# lines(log(n_seq)[index_set],log(rowMeans((abs(n_record_theta_SGP-L_2_minimizer)))),type='p',col='blue',
#       pch=1,cex=0.8)
lines(log(n_seq)[index_set],log((n_seq^{-3/(2*3+1)}/38))[index_set],type='l',lty=1,cex=0.8)
#lines(log(n_seq)[index_set],log((n_seq^{-1/(2)}/22))[index_set],type='l',lty=1,cex=0.8)

legend("bottomleft", 
       legend = c("GaSP", "S-GaSP"), 
       col = c('red','blue'),
       pch=c(2,1),cex=0.8)
dev.off()



