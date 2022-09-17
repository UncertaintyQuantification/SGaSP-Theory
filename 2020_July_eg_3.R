

library(RobustCalibration)
library(RobustGaSP)
library(nloptr)
library(lhs)
library(DiceKriging)
library(plot3D)




reality<- function(xx)
{
  x1 <- xx[,1]
  x2 <- xx[,2]
  
  # 2/3*exp(x1+0.2)-x2*sin(0.4)+0.4+10*exp(-x1)*(x1-1/2)*(x2^2-x2+1/6)
  
  ##CHENG & SANDU (2010) FUNCTION
  #https://www.sfu.ca/~ssurjano/chsan10.html
  #cos(x1+x2)*exp(x1*x2)
  
  (sin(0.2*pi*x1)+sin(2*pi*x1))*x2+1
}

cm_model <- function(x, theta)
{
  x1 <- x[,1]
  x2 <- x[,2]
  
  y=sin(theta*x1)*x2
  #y=5*(x-theta)^2
  
  return(y)
}



cm_model_mean <- function(x, theta){
  #y=sin(theta[1]*x-theta[2])+theta[3]*x^2+theta[4]
  
  x1 <- x[,1]
  x2 <- x[,2]
  
  y=sin(theta[1]*x1)*x2+theta[2]
  return(y)
}
LS_calibration<-function(theta,X_here){
  output_tilde=output-cm_model(input,theta)
  
  LS=(t(output_tilde) - t(output_tilde)%*%X_here%*%solve(t(X_here)%*%X_here)%*%t(X_here))%*%output_tilde
  LS
  
}

LS_calibration_2<-function(theta){
  
 mean((cm_model_mean(input,theta)-output)^2)
}


##this one give closed form expression for the mean
L2_calibration<-function(theta,pred_testing_output,X_testing_here){
  
  pred_testing_output_tilde=pred_testing_output-cm_model(testing_input,theta)
  
  L2_loss=(t(pred_testing_output_tilde) - t(pred_testing_output_tilde)%*%X_testing_here%*%solve(t(X_testing_here)%*%X_testing_here)%*%t(X_testing_here))%*%pred_testing_output_tilde
  L2_loss
}


L2_calibration_2<-function(theta,pred_testing_output){
  mean((cm_model_mean(testing_input,theta)-pred_testing_output)^2)
}



neg_profile_lik_all<-function(param,have_cm_model=T,have_trend=T,output_here=output){
  
  if(have_cm_model==T){
    theta=(param[1:p_theta])
  }
  beta=exp(param[(p_theta+1):(p_theta+p_x)])
  
  eta=exp(param[p_theta+p_x+1])
  
  R=separable_kernel(R0,beta,kernel_type='matern_5_2',alpha=rep(2,p_x))
  
  #separable_kernel(R0, beta, kernel_type, alpha)
  
  R_tilde=R+eta*diag(n)
  
  L=t(chol(R_tilde))
  
  output_tilde=output_here
  
  if(have_cm_model==T){
    output_tilde=output_here-cm_model(input,theta)
  }
  
  ###if there is mean parameter
  if(have_trend==T){
    R_inv_X=(backsolve(t(L),forwardsolve(L,X)))
    L_X=t(chol(t(X)%*%R_inv_X))
    X_t_R_inv_y=t(R_inv_X)%*%output_tilde
    
    theta_hat=(backsolve(t(L_X),forwardsolve(L_X,X_t_R_inv_y)))
    output_tilde=output_tilde-X%*%theta_hat
  }
  ###
  
  R_inv_output_tilde=(backsolve(t(L),forwardsolve(L,output_tilde)))
  
  S_2=t(output_tilde)%*%R_inv_output_tilde
  
  
  -(-sum(log(diag(L)))-length(output_here)/2*log(S_2))
  
}



M=100



RMSE_all_record=matrix(0,M,8) 
RMSE_fm_record=matrix(0,M,8) 

##
avg_interval_record=matrix(0,M,8) 
avg_prop_record=matrix(0,M,8) 

q=1
p_theta=1
p_x=1


record_par_gasp=matrix(0,M,p_theta+1)
record_par_sgasp=matrix(0,M,p_theta+1)



record_par_l2=matrix(0,M,p_theta+1)
record_par_ls=matrix(0,M,p_theta+1)
record_par_no_discrepancy=matrix(0,M,p_theta+1)

record_par_gasp_mle=matrix(0,M,p_theta+1)
record_par_sgasp_mle=matrix(0,M,p_theta+1)
record_par_no_discrepancy_mle=matrix(0,M,p_theta+1)


p_x=2
p_theta=1

n=50

time_record=system.time(
  for(i_M in 1:M){
    print(i_M)
    set.seed(i_M)
    input=maximinLHS(n=n,k=p_x)

    have_cm_model=T
    have_trend=T
    n=dim(input)[1]
    output=reality(input)
    
    output=reality(input)+rnorm(n,mean=0,sd=0.1)
    #X=cbind(rep(1,n),input)
    X=as.matrix(rep(1,n))
    
    q=dim(X)[2]
    p_theta=1
    
    
    n_sub_testing=50
    n_testing=n_sub_testing^2
    testing_input1=as.vector(seq(0,1,1/(n_sub_testing-1)))
    testing_input2=as.vector(seq(0,1,1/(n_sub_testing-1)))
    
    
    testing_input=cbind( rep(seq(0,1,1/(n_sub_testing-1)),n_sub_testing),
                         as.vector(t(matrix(seq(0,1,1/(n_sub_testing-1)),n_sub_testing,n_sub_testing))))
    
    reality_testing=reality(testing_input)
    testing_output=reality_testing+rnorm(n_testing,mean=0,sd=0.1)
    #as.matrix(seq(0,1,1/(n_testing-1)))
    X_testing=as.matrix(rep(1,n_testing))
    
    
    
    ###gasp, posterior sample
    
    model_gasp=rcalibration(design=input, observations=output, simul_type=1,
                             have_trend=T,X=X,
                             math_model=cm_model,theta_range=matrix(c(0,10),1,2)
                             ,S=12000,S_0=2000,discrepancy_type='GaSP')
    
    model_gasp_pred=predict(model_gasp,testing_input,X_testing=X_testing,math_model=cm_model,
                             interval_est=c(0.025,0.975),interval_data=F,
                             n_thinning =10 )
    
    
    RMSE_all_record[i_M,1]=sqrt(mean( (model_gasp_pred@mean- reality_testing)^2))
    
    RMSE_fm_record[i_M,1]=sqrt(mean( (model_gasp_pred@math_model_mean- reality_testing)^2))
    
    avg_interval_record[i_M,1]=mean(model_gasp_pred@interval[,2]-model_gasp_pred@interval[,1])
    
    avg_prop_record[i_M,1]=length(which( (reality_testing> model_gasp_pred@interval[,1]) & (reality_testing<model_gasp_pred@interval[,2]) ))/n_testing
    
    record_par_gasp[i_M,]=colMeans(model_gasp@post_sample)[c(1,6)]

    
    
    ###sgasp
    model_sgasp=rcalibration(design=input, observations=output, simul_type=1,
                             have_trend=T,X=X,
                             math_model=cm_model,theta_range=matrix(c(0,10),1,2)
                             ,S=12000,S_0=2000,discrepancy_type='S-GaSP')

    ##interval for the mean
    model_sgasp_pred=predict(model_sgasp,testing_input,X_testing=X_testing,math_model=cm_model,
                             interval_est=c(0.025,0.975),interval_data=F,
                             n_thinning =10 )
    
    RMSE_all_record[i_M,2]=sqrt(mean( (model_sgasp_pred@mean- reality_testing)^2))
    RMSE_fm_record[i_M,2]=sqrt(mean( (model_sgasp_pred@math_model_mean- reality_testing)^2))
    
    avg_interval_record[i_M,2]=mean(model_sgasp_pred@interval[,2]-model_sgasp_pred@interval[,1])
    
    avg_prop_record[i_M,2]=length(which( (reality_testing> model_sgasp_pred@interval[,1]) & (reality_testing<model_sgasp_pred@interval[,2]) ))/n_testing
    
    
    record_par_sgasp[i_M,]=colMeans(model_sgasp@post_sample)[c(1,6)]
    #median(model_sgasp@post_sample[,1])
    #image2D(matrix(reality_testing,n_sub_testing,n_sub_testing),testing_input1,testing_input2)
    #image2D(matrix(model_sgasp_pred@mean,n_sub_testing,n_sub_testing),testing_input1,testing_input2)
    
    
    #model_sgasp_pred@mean
    
    ###two steps, L2, 
    m_l2=rgasp(input,output,nugget.est=T,method='mle')
    pred_l2_all=predict(m_l2,testing_input,interval_data=F)
    pred_l2=pred_l2_all$mean
    
    avg_interval_record[i_M,3]=mean(pred_l2_all$upper95-pred_l2_all$lower95)
    
    avg_prop_record[i_M,3]=length(which( (reality_testing> pred_l2_all$lower95) & (reality_testing<pred_l2_all$upper95) ))/n_testing
    
    

    #this is profile likelihood by taking hat theta_2 into account
    
    m_l2_est=try(optim(c(1), L2_calibration,pred_testing_output=pred_l2,X_testing_here=X_testing,method="Brent",lower=0,upper=10 ),silent = T)
    

    theta_hat_l2=solve(t(X_testing)%*%X_testing)%*%t(X_testing)%*%(pred_l2-cm_model(testing_input,m_l2_est$par))
    
    RMSE_fm_record[i_M,3]=sqrt(mean( (cm_model(testing_input,m_l2_est$par)+X_testing%*%theta_hat_l2-reality_testing)^2))
    RMSE_all_record[i_M,3]=sqrt(mean( (pred_l2-reality_testing)^2))
    
    record_par_l2[i_M,]=c(m_l2_est$par,theta_hat_l2)
    

    ###ls
    m_ls_record=try(optim(c(5), LS_calibration, X_here=X, method="Brent",lower=0,upper=10),silent = T)
    
    theta_hat_ls= solve(t(X)%*%X)%*%t(X)%*%(output-cm_model(input,m_ls_record$par))
    output_tilde_ls=output-cm_model(input,m_ls_record$par)-X%*%theta_hat_ls
    
    m_ls=rgasp(input,output_tilde_ls,nugget.est=T,method='mle')
    pred_ls_all=predict(m_ls,testing_input,interval_data=F)
    
    pred_ls_mean=predict(m_ls,testing_input)$mean+cm_model(testing_input,m_ls_record$par)+X_testing%*%theta_hat_ls
    
    pred_ls_all_upper95=pred_ls_all$upper95+cm_model(testing_input,m_ls_record$par)+X_testing%*%theta_hat_ls
    pred_ls_all_lower95=pred_ls_all$lower95+cm_model(testing_input,m_ls_record$par)+X_testing%*%theta_hat_ls
                                                                                      
    avg_interval_record[i_M,4]=mean(pred_ls_all_upper95-pred_ls_all_lower95)
    
    avg_prop_record[i_M,4]=length(which( (reality_testing> (pred_ls_all_lower95) ) & (reality_testing<pred_ls_all_upper95) ))/n_testing
    
    RMSE_fm_record[i_M,4]=sqrt(mean((cm_model(testing_input,m_ls_record$par)+X_testing%*%theta_hat_ls-reality_testing)^2))
    
    RMSE_all_record[i_M,4]=sqrt(mean((pred_ls_mean-reality_testing)^2))
    
    record_par_ls[i_M,]=c(m_ls_record$par,theta_hat_ls)
    
    m_no_discrepancy=rcalibration(input,output,p_theta=p_theta,simul_type=1,math_model=cm_model,have_trend=T,X=X,
                                  theta_range=matrix(c(0,10),1,2),S=12000,S_0=2000,discrepancy_type = 'no-discrepancy')
    #plot(m_no_discrepancy@post_sample[,1])
    #den_est_no_dis=(kde2d(m_no_discrepancy@post_sample[,1],m_no_discrepancy@post_sample[,2]))
    #contour(den_est_no_dis, xlab = "theta1",
    #        ylab = "theta2" ,xlim=theta_range[1,],ylim=theta_range[2,])
    
    
    
    
    m_no_discrepancy_pred=predict(m_no_discrepancy,testing_input,math_model=cm_model,n_thinning=10,
                                  interval_est=c(0.025,0.975),interval_data=F,X_testing=X_testing)
    
    RMSE_fm_record[i_M,5]=sqrt(mean( (m_no_discrepancy_pred@math_model_mean- reality_testing)^2))
    
    avg_interval_record[i_M,5]=mean(m_no_discrepancy_pred@interval[,2]-m_no_discrepancy_pred@interval[,1])
    
    avg_prop_record[i_M,5]=length(which( (reality_testing> m_no_discrepancy_pred@interval[,1]) & (reality_testing<m_no_discrepancy_pred@interval[,2]) ))/n_testing
    
    record_par_no_discrepancy[i_M,]=colMeans(m_no_discrepancy@post_sample)[c(1,3)]
    
    

    ###mle, gasp
    model_gasp_mle=rcalibration(design=input, observations=output, simul_type=1,
                                have_trend=T,X=X,
                                math_model=cm_model,theta_range=matrix(c(0,10),1,2),discrepancy_type='GaSP',method='mle',num_initial_values=5)
    
    model_gasp_pred_mle=predict(model_gasp_mle,testing_input,X_testing=X_testing,math_model=cm_model,
                                interval_est=c(0.025,0.975),interval_data=F)
    
    
    RMSE_all_record[i_M,6]=sqrt(mean( (model_gasp_pred_mle@mean- reality_testing)^2))
    
    RMSE_fm_record[i_M,6]=sqrt(mean( (model_gasp_pred_mle@math_model_mean- reality_testing)^2))
    
    avg_interval_record[i_M,6]=mean(model_gasp_pred_mle@interval[,2]-model_gasp_pred_mle@interval[,1])
    
    avg_prop_record[i_M,6]=length(which( (reality_testing> model_gasp_pred_mle@interval[,1]) & (reality_testing<model_gasp_pred_mle@interval[,2]) ))/n_testing
    
    record_par_gasp_mle[i_M,]=(model_gasp_mle@param_est)[c(1,6)] ##the six is the mean
    
    
    ##mle, sgasp
    model_sgasp_mle=rcalibration(design=input, observations=output, simul_type=1,
                                have_trend=T,X=X,
                                math_model=cm_model,theta_range=matrix(c(0,10),1,2),discrepancy_type='S-GaSP',method='mle',num_initial_values=5)
    
    model_sgasp_pred_mle=predict(model_sgasp_mle,testing_input,X_testing=X_testing,math_model=cm_model,
                                interval_est=c(0.025,0.975),interval_data=F )
    
    
    RMSE_all_record[i_M,7]=sqrt(mean( (model_sgasp_pred_mle@mean- reality_testing)^2))
    
    RMSE_fm_record[i_M,7]=sqrt(mean( (model_sgasp_pred_mle@math_model_mean- reality_testing)^2))
    
    avg_interval_record[i_M,7]=mean(model_sgasp_pred_mle@interval[,2]-model_sgasp_pred_mle@interval[,1])
    
    avg_prop_record[i_M,7]=length(which( (reality_testing> model_sgasp_pred_mle@interval[,1]) & (reality_testing<model_sgasp_pred_mle@interval[,2]) ))/n_testing
    
    record_par_sgasp_mle[i_M,]=(model_sgasp_mle@param_est)[c(1,6)] ##the six is the mean
    
    
    ###mle, robustcalibration
    model_no_discrepancy_mle=rcalibration(design=input, observations=output, simul_type=1,
                                have_trend=T,X=X,
                                math_model=cm_model,theta_range=matrix(c(0,10),1,2),discrepancy_type='no-discrepancy',method='mle',num_initial_values=5)
    
    model_no_discrepancy_pred_mle=predict(model_no_discrepancy_mle,testing_input,X_testing=X_testing,math_model=cm_model,
                                interval_est=c(0.025,0.975),interval_data=F)
    
    
    #RMSE_all_record[i_M,]=sqrt(mean( (model_gasp_pred_mle@mean- reality_testing)^2))
    
    RMSE_fm_record[i_M,8]=sqrt(mean( (model_no_discrepancy_pred_mle@math_model_mean- reality_testing)^2))
    
    #avg_interval_record[i_M,8]=mean(model_no_discrepancy_pred_mle@interval[,2]-model_no_discrepancy_pred_mle@interval[,1])
    
    #avg_prop_record[i_M,8]=length(which( (reality_testing> model_no_discrepancy_pred_mle@interval[,1]) & (reality_testing<model_no_discrepancy_pred_mle@interval[,2]) ))/n_testing
    
    record_par_no_discrepancy_mle[i_M,]=(model_no_discrepancy_mle@param_est)[c(1,3)] ##the six is the mean
    
    
    
     print(RMSE_fm_record[i_M,])
     print(RMSE_all_record[i_M,])
    # if(i_M>1){
    #   print(colMeans(avg_prop_record[i_M,]))
    # }
    # 
  }
  

  
)

# pdf('reality_L2_SGaSP_eg3.pdf',width=10,height=3.5)
# par(mfrow=c(1,3),mai=c(.7,.8,.8,.28))
#  image2D( matrix(abs(reality_testing),50,50),zlim=c(-0.35,2.3),,xlab=expression(x[1]),ylab=expression(x[2]))
# # image2D( matrix(abs(pred_l2),50,50),zlim=c(-0.35,2.3))
# # image2D( matrix(abs(model_sgasp_pred@mean),50,50),zlim=c(-0.35,2.3))
#  image2D( matrix((pred_l2-reality_testing),50,50),zlim=c(-0.33,0.33),xlab=expression(x[1]),ylab=expression(x[2]))
# ### image2D( matrix((pred_ls_mean-reality_testing),50,50),zlim=c(-0.33,0.33))
#  image2D( matrix((model_sgasp_pred@mean-reality_testing),50,50),zlim=c(-0.33,0.33),xlab=expression(x[1]),ylab=expression(x[2]))
# dev.off()



colMeans(RMSE_fm_record)
colMeans(RMSE_all_record)
colMeans(avg_interval_record)
colMeans(avg_prop_record)

pdf("boxplot_eg_3_theta_1_n_50.pdf",width=6,height=7)
par(mgp=c(3, 1.8, 0),  las=1)
theta_1=cbind(record_par_gasp[,1],record_par_gasp_mle[,1],record_par_sgasp[,1],record_par_sgasp_mle[,1],record_par_l2[,1],record_par_ls[,1])
boxplot(theta_1,names=c('GaSP \n  PS','GaSP\n MLE','S-GaSP\n   PS','S-GaSP \n MLE',expression(L[2]), 'LS'),
        ylab='',main=expression(theta[1])) 
abline(a=6.494727,b=0,lty=1)
dev.off()

pdf("boxplot_eg_3_theta_2_n_50.pdf",width=6,height=7)
par(mgp=c(3, 1.8, 0),  las=1)
theta_2=cbind(record_par_gasp[,2],record_par_gasp_mle[,2],record_par_sgasp[,2],record_par_sgasp_mle[,2],record_par_l2[,2],record_par_ls[,2])
boxplot(theta_2,names=c('GaSP \n  PS','GaSP\n MLE','S-GaSP\n   PS','S-GaSP \n MLE',expression(L[2]), 'LS'),
        ylab='',main=expression(theta[2])) 
abline(a=1.149148,b=0,lty=1)
dev.off()

