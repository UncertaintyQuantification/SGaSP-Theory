##rewrite this code by the  code from FastGaSP package and a new function for linear coefficient
#library(FastGaSP)

library(Rcpp)
library(RcppEigen)
library(RobustGaSP)
sourceCpp(file='src/functions_June_30_2020.cpp') 


y_R<-function(x){
  j_seq=seq(1,100,1)
  record_y_R=0
  for(i_j in 1:100){
    record_y_R=record_y_R+1*j_seq[i_j]^{-1}*sin(j_seq[i_j]*5)*cos(pi*(j_seq[i_j]-0.5)*x*5)+1*j_seq[i_j]^{-2*3}*sin(j_seq[i_j]*5)*cos(pi*(j_seq[i_j]-0.5)*x*5)
    
  }
  record_y_R
}


f_M<-function(x,theta){
  j_seq=seq(1,100,1)
  record_f_M=0
  for(i_j in 1:100){
    #record_y_R=record_y_R+2*j_seq[i_j]^{-2*3}*sin(j_seq[i_j])*cos(pi*(j_seq[i_j]-0.5)*x)
    record_f_M=record_f_M+theta*j_seq[i_j]^{-1}*sin(j_seq[i_j]*5)*cos(pi*(j_seq[i_j]-0.5)*x*5)
    #record_y_R=record_y_R+2*j_seq[i_j]^{-2*2}*sin(j_seq[i_j]*20)*cos(pi*(j_seq[i_j]-0.5)*x*20)
    
  }
  record_f_M
}

#mean((output-f_M(input,1))^2)
#mean((output-f_M(input,4))^2)

plot(y_R(seq(0,1,0.001)))

lines((f_M(seq(0,1,0.001),3)),type='l')

plot(y_R(seq(0,1,0.001))-(f_M(seq(0,1,0.001),1)),type='l')

mean( (y_R(seq(0,1,0.001))-(f_M(seq(0,1,0.001),1)))^2)
mean( (y_R(seq(0,1,0.001))-(f_M(seq(0,1,0.001),2)))^2)

pred_mean_FFBS_X_GP<-function(sigma_2,beta,eta, delta_x,index_obs,output){
  
  param=log(c(beta,eta))
  have_noise=T
  theta_hat=Get_linear_coeff_1_dim_GP(param, have_noise, delta_x, output,
                                      kernel_type,  X)

  return_list=as.list(1:2)
  
  return_list[[1]]=theta_hat
  return_list[[2]]=X_testing*theta_hat+Kalman_smoother(param, have_noise, index_obs, 
                                                       delta_x_all,  output-X*theta_hat, 1.0,
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
                                                                   delta_x_all,  output-X*theta_hat, 1.0,
                                                                   kernel_type,only_mean)[[1]]  
  return_list
}


log_lik_zero_mean<-function(param,have_noise,delta_x, output_here, kernel_type){
  param=as.vector(param)
  
  log_det_S2=Get_log_det_S2(param,have_noise,delta_x,
                            output_here, kernel_type);
  ##log likelihood
  -log_det_S2[[1]]/2-(num_obs)/2*log(log_det_S2[[2]])
}


#n_sim=50
#n_experiment=100
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
have_noise=T

#n_record_theta=matrix(0,n_sim,2)
#n_record_RMSE=matrix(0,n_sim,2)

num_obs_all=3*10^4



n_record_theta_LS=matrix(0,n_sim,n_experiment)
n_record_RMSE_LS=matrix(0,n_sim,n_experiment)

param_L2_record=array(0,c(n_sim,n_experiment,2))

n_record_theta_L2=matrix(0,n_sim,n_experiment)
n_record_RMSE_L2=matrix(0,n_sim,n_experiment)

param_LS_record=array(0,c(n_sim,n_experiment,2))




param_GP_record=array(0,c(n_sim,n_experiment,2))
param_SGP_record=array(0,c(n_sim,n_experiment,2))

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
    X=rep(0,num_obs)
    j_seq=seq(1,100,1)
    for(i_j in 1:100){

      X=X+j_seq[i_j]^{-1}*sin(j_seq[i_j]*5)*cos(pi*(j_seq[i_j]-0.5)*input*5)

    }
    
    X_testing=rep(0,num_obs_all)
    for(i_j in 1:100){
      X_testing=X_testing+j_seq[i_j]^{-1}*sin(j_seq[i_j]*5)*cos(pi*(j_seq[i_j]-0.5)*input_all*5)
    }
    
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
      sigma_2=1 ##the result won't change for different sigma_2

      model_GP=optim(c(0,-5), log_lik_GP_with_1_dim_X, have_noise=have_noise,delta_x=delta_x,output=output,
                         kernel_type=kernel_type,X=X,control = list(fnscale=-1,maxit=100))
      

      param_GP= model_GP$par
      pred_mean_X_all_GP=pred_mean_FFBS_X_GP(sigma_2,exp(param_GP[1]),exp(param_GP[2]), delta_x,index_obs,output)
      

      n_record_theta_GP[i_n,i_exp]=pred_mean_X_all_GP[[1]]  ##theta
      n_record_RMSE_GP[i_n,i_exp]=sqrt(mean((pred_mean_X_all_GP[[2]]-testing_output)^2))
      

      
      model_SGP=optim(c(0,-5), log_lik_SGP_with_1_dim_X, have_noise=have_noise,delta_x=delta_x,output=output,
                     kernel_type=kernel_type,X=X,control = list(fnscale=-1,maxit=100))
      
      
      param_SGP= model_SGP$par
      
      lambda_z_est=sqrt(num_obs*exp(param_SGP[1])/exp(param_SGP[2]));
      
      lambda_est=exp(param_SGP[2])/num_obs
      shrinkage=1+lambda_est*lambda_z_est
      
      pred_mean_X_all_SGP=pred_mean_FFBS_X_SGP(sigma_2,exp(param_SGP[1]),exp(param_SGP[2]),lambda_z_est, delta_x,index_obs,output)
      
      n_record_theta_SGP[i_n,i_exp]=pred_mean_X_all_SGP[[1]]  ##theta
      
      n_record_RMSE_SGP[i_n,i_exp]=sqrt(mean((pred_mean_X_all_SGP[[2]]-testing_output)^2))
      
      param_GP_record[i_n,i_exp,]=exp(model_GP$par)
      param_SGP_record[i_n,i_exp,]=exp(model_SGP$par)
      
      
      ###LS, L2
      theta_hat_LS=1/(sum(X^2))*(t(X)%*%output)
      output_tilde_LS=output-X%*%theta_hat_LS
      
      model_LS=optim(c(0,-5), log_lik_zero_mean, have_noise=have_noise,delta_x=delta_x,output=output,
                     kernel_type=kernel_type,control = list(fnscale=-1,maxit=100))
      
      param_LS=model_LS$par
      
      
      pred_LS=X_testing%*%theta_hat_LS+Kalman_smoother(param_LS, have_noise, index_obs, 
                                                       delta_x_all,  output_tilde_LS, 1.0,
                                                       kernel_type,only_mean)[[1]]  
      
      
      n_record_theta_LS[i_n,i_exp]=theta_hat_LS  ##theta
      
      n_record_RMSE_LS[i_n,i_exp]=sqrt(mean((pred_LS-testing_output)^2))
      
      param_LS_record[i_n,i_exp,]=exp(param_LS)
      
      
      ###L2
      model_L2=optim(c(0,-5), log_lik_zero_mean, have_noise=have_noise,delta_x=delta_x,output=output,
                     kernel_type=kernel_type,control = list(fnscale=-1,maxit=100))
      
      pred_L2=Kalman_smoother(param_LS, have_noise, index_obs, 
                              delta_x_all,  output, 1.0,
                              kernel_type,only_mean)[[1]]  
      theta_L2=1/(sum(X_testing^2))*(t(X_testing)%*%pred_L2)
      
      n_record_theta_L2[i_n,i_exp]=theta_L2  ##theta
      
      n_record_RMSE_L2[i_n,i_exp]=sqrt(mean((pred_L2-testing_output)^2))
      
      param_LS_record[i_n,i_exp,]=exp(model_L2$par)
      
      
    }
    print(c( mean(n_record_RMSE_GP[i_n,]),mean(n_record_RMSE_SGP[i_n,]) ))
    
    print(c( mean(n_record_RMSE_LS[i_n,]),mean(n_record_RMSE_L2[i_n,]) ))
    

  }
)



n_record_RMSE_GP_mean=matrix(0,n_sim,n_experiment )
n_record_RMSE_SGP_mean=matrix(0,n_sim,n_experiment )
n_record_RMSE_LS_mean=matrix(0,n_sim,n_experiment )
n_record_RMSE_L2_mean=matrix(0,n_sim,n_experiment )

for(i_sim in 1:n_sim){
  for(i_experiment in 1:n_experiment){
    n_record_RMSE_GP_mean[i_sim,i_experiment]=sqrt( mean((testing_output-n_record_theta_GP[i_sim,i_experiment])^2))
    n_record_RMSE_SGP_mean[i_sim,i_experiment]=sqrt( mean((testing_output-n_record_theta_SGP[i_sim,i_experiment])^2))
    n_record_RMSE_LS_mean[i_sim,i_experiment]=sqrt( mean((testing_output-n_record_theta_LS[i_sim,i_experiment])^2))
    n_record_RMSE_L2_mean[i_sim,i_experiment]=sqrt( mean((testing_output-n_record_theta_L2[i_sim,i_experiment])^2))
    
  }
}

n_seq=floor(exp(5+(1:n_sim)/10))


L_2_minimizer= as.numeric(solve(t(X_testing)%*% X_testing)%*%t(X_testing)%*%(testing_output))

index_set=1:50

library(ggplot2)


record_RMSE_est_all=data.frame(method = factor(rep(c("GaSP", "S-GaSP","L2","LS"), each=n_sim)), 
                                  RMSE = log((c(as.vector(rowMeans(n_record_RMSE_GP) ),as.vector(rowMeans(n_record_RMSE_SGP)),
                                           as.vector(rowMeans(n_record_RMSE_L2)),as.vector(rowMeans(n_record_RMSE_LS))))) )

pdf("pred_GP_SGP_L2_LS_eg2.pdf",heigh=3,width=4.5)

ggplot() + 
  geom_point(data = record_RMSE_est_all, aes( x=rep(log(n_seq),4), y = RMSE, shape=method,color = method),size=1) +
  labs(x ='log(n)',y=expression(log(AvgRMSE[f^M+delta]) ))+
  scale_color_manual(values = c("GaSP" = "red", "S-GaSP" = "blue","L2"="green","LS"="brown"))+
  scale_shape_manual(values=c("GaSP" =17, "S-GaSP" = 16,"L2"=15,"LS"=18))
  
dev.off()




record_theta_est_GP_SGP=data.frame(method = factor(rep(c("GaSP", "S-GaSP"), each=n_sim*n_experiment)), 
                            theta = c(as.vector(n_record_theta_GP),as.vector(n_record_theta_SGP)))


ggplot(record_theta_est_GP_SGP, aes(x=theta , fill=method)) +
  geom_histogram(binwidth=.002, alpha=.5, position="identity")+labs(x = expression(theta))




pdf("hist_GP_SGP_eg2.pdf",heigh=3,width=4.5)

ggplot(record_theta_est_GP_SGP, aes(x=theta , color=method)) +
  geom_histogram(binwidth=.002, alpha=.5, position="identity",fill='white')+
  labs(x =expression(theta))+ scale_colour_manual(name = "method", values = c("red", "blue"), 
                                                     labels = c('GaSP', 'S-GaSP')) 

dev.off()
record_theta_est_L2_LS=data.frame(method = factor(rep(c("L2", "LS"), each=n_sim*n_experiment)), 
                                   theta = c(as.vector(n_record_theta_L2),as.vector(n_record_theta_LS)))


ggplot(record_theta_est_L2_LS, aes(x=theta , fill=method)) +
  geom_histogram(binwidth=.002, alpha=.5, position="identity")+labs(x = expression(theta))

pdf("hist_L2_LS_eg2.pdf",heigh=3,width=4.5)

ggplot(record_theta_est_L2_LS, aes(x=theta , color=method)) +
  geom_histogram(binwidth=.002, alpha=.5, position="identity",fill='white')+
  labs(x =expression(theta))+scale_colour_manual(name = "method", values = c("green", "brown"), 
                        labels = c(expression(L[2]), 'LS') ) 
dev.off()


