library(deSolve)
library(RobustCalibration)
library(RobustGaSP)
library(MASS)
library(plot3D)



Box_model <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dM1=-10^{parameters[1]-3}*M1
    dM2=10^{parameters[1]-3}*M1-10^{parameters[2]-3}*M2
    list(c(dM1, dM2))
  })
}

Box_model_solved<-function(input, theta){
  init <- c(M1 = 100 , M2 = 0)
  out <- ode(y = init, times = c(0,input), func = Box_model, parms = theta)
  return(out[-1,3])
}



LS_calibration<-function(theta,input_here,output_here){
  mean((Box_model_solved(input_here,theta)-output_here)^2)
}

L2_calibration<-function(theta,testing_input_here,pred_testing_output){
  mean((Box_model_solved(testing_input_here,theta)-pred_testing_output)^2)
}


output_all=matrix(c(19.2,14,14.4,24,42.3,30.8,42.1,40.5,40.7,46.4,27.1,22.3))
input_all=matrix(c(10,10,20,20,40,40,80,80,160,160,320,320))

n=6
k=2
output_mat_all=t(matrix(output_all,k,n))
input_mat_all=matrix(c(10,20,40,80,160,320))


input_plot=as.matrix(input_all)
output_plot=as.matrix(output_all)
theta_test=c(1,.8)
plot(input_plot,Box_model_solved(input_plot,theta_test),type='l',ylim=c(0,50))
lines(input_plot,output_plot,type='p',pch=20)

set.seed(1)
##fit 
output=as.matrix(output_mat_all)
input=as.matrix(input_mat_all)

input_seq=as.matrix(seq(1,350,1))

theta_range=matrix(c(0.5,0.5,1.5,1.5),2,2)

m_no_discrepancy=rcalibration(input,output,simul_type=1,math_model=Box_model_solved,
                              sd_proposal=c(0.25,0.25,1,1),
                              theta_range=theta_range,S=12000,S_0=2000,discrepancy_type = 'no-discrepancy')

m_no_discrepancy_pred=predict(m_no_discrepancy,input_seq,math_model=Box_model_solved,n_thinning=5,
                              interval_est=c(0.025, 0.975),interval_data=T)

m_gasp=rcalibration(input,output,simul_type=1,math_model=Box_model_solved,
                    sd_proposal=c(0.25,0.25,1,1),
                    theta_range=theta_range,S=12000,S_0=2000,discrepancy_type = 'GaSP')


m_gasp_pred=predict(m_gasp,input_seq,math_model=Box_model_solved,n_thinning=5,
                    interval_est=c(0.025, 0.975),interval_data=T)

m_sgasp=rcalibration(input,output,simul_type=1,math_model=Box_model_solved,
                     sd_proposal=c(0.25,0.25,1,1),
                     theta_range=theta_range,S=12000,S_0=2000)


model_sgasp_pred=predict(m_sgasp,input_seq,math_model=Box_model_solved,n_thinning=5,
                         interval_est=c(0.025, 0.975),interval_data=T)


###two steps
output_vec=as.vector(t(output))
input_vec=as.vector(t(matrix(input,length(input),2)))
##l2
m_l2=rgasp(input_vec,output_vec,nugget.est=T)
pred_l2=predict(m_l2,input_seq,interval_data=T)



m_l2_record=try(optim(c(1,1), L2_calibration, testing_input_here=input_seq, pred_testing_output=pred_l2$mean, 
                      method="L-BFGS-B",lower=theta_range[,1],upper=theta_range[,2]),silent = T)

pred_l2_cm=Box_model_solved(input_seq,m_l2_record$par)

##ls
m_ls_record=try(optim(c(1,1), LS_calibration,input_here=input_vec,output_here=output_vec,
                      method="L-BFGS-B",lower=theta_range[,1],upper=theta_range[,2]),silent = T)


output_tilde_ls=output_vec-Box_model_solved(input_vec,m_ls_record$par)

m_ls=rgasp(input_vec,output_tilde_ls,nugget.est=T)
pred_ls_cm_mean=Box_model_solved(input_seq,m_ls_record$par)
m_ls_pred=predict(m_ls,input_seq)
pred_ls_mean=m_ls_pred$mean+pred_ls_cm_mean

###mle 
m_gasp_mle=rcalibration(design=input, observations=output,simul_type=1,math_model=Box_model_solved,
                        theta_range=theta_range,discrepancy_type = 'GaSP',method='mle',num_initial_values=5)

m_sgasp_mle=rcalibration(design=input, observations=output,simul_type=1,math_model=Box_model_solved,
                         theta_range=theta_range,discrepancy_type = 'S-GaSP',method='mle',num_initial_values=5)
m_no_discrepancy_mle=rcalibration(design=input, observations=output,simul_type=1,math_model=Box_model_solved,
                                  theta_range=theta_range,discrepancy_type = 'no-discrepancy',method='mle',num_initial_values=5)



m_gasp_mle_predict=predict(m_gasp_mle,input_seq,math_model=Box_model_solved,
                           interval_est=c(0.025, 0.975),interval_data=T)  ##no mean, maybe set mean to be math model mean?

m_sgasp_mle_predict=predict(m_sgasp_mle,input_seq,math_model=Box_model_solved,
                            interval_est=c(0.025, 0.975),interval_data=T)  ##no mean, maybe set mean to be math model mean?

m_no_discrepancy_mle_predict=predict(m_no_discrepancy_mle,input_seq,math_model=Box_model_solved,
                                     interval_est=c(0.025, 0.975),interval_data=T)  ##no mean, maybe set mean to be math model mean?



pdf('plot_cm_pred_eg_4_box.pdf',height=4.5,width=6)

plot(input_seq,m_gasp_pred@mean,type='l',col='red',ylim=c(0,60),
     xlab='t',ylab=expression(y[2]),main='Calibrated computer model')
#
polygon( c(input_seq,rev(input_seq)),c(m_gasp_pred@interval[,1],
                                       rev(m_gasp_pred@interval[,2])),col =  "grey80", border = FALSE)
lines(input_seq,m_no_discrepancy_pred@math_model_mean,type='l',col='green',lty=1)
lines(input_seq,model_sgasp_pred@math_model_mean,type='l',col='blue',lty=1)
lines(input_seq,m_gasp_pred@math_model_mean,type='l',col='red',lty=1)

lines(input_seq,m_no_discrepancy_mle_predict@math_model_mean,type='l',col='green',lty=2)
lines(input_seq,m_sgasp_mle_predict@math_model_mean,type='l',col='blue',lty=2)
lines(input_seq,m_gasp_mle_predict@math_model_mean,type='l',col='red',lty=2)

lines(input_seq,pred_l2_cm,col='brown',lty=2)
lines(input_seq,pred_ls_cm_mean,col='orange',lty=2)

lines(input_all,output_all,type='p',col='black',pch=20)
legend("bottom", legend=c('N-D, PS','S-GaSP, PS','GaSP, PS',expression(L[2]),
                            'N-D, MLE','S-GaSP, MLE','GaSP, MLE','LS'),
       lty=c(1,1,1,2,2,2,2,2),col=c('green','blue','red','brown','green','blue','red','orange'),ncol=2, x.intersp=1,cex=.8)

dev.off()


pdf('plot_pred_eg_4_box.pdf',height=4.5,width=6)

plot(input_seq,m_gasp_pred@mean,type='l',col='red',ylim=c(0,60),xlab='t',ylab=expression(y[2]),
     main='Calibrated computer model+discrepancy')
#
polygon( c(input_seq,rev(input_seq)),c(model_sgasp_pred@interval[,1],
                                       rev(model_sgasp_pred@interval[,2])),col =  "grey80", border = FALSE)
lines(input_seq,m_gasp_pred@mean,type='l',col='red')
lines(input_seq,model_sgasp_pred@mean,type='l',col='blue',lty=1)
lines(input_seq,m_gasp_mle_predict@mean,type='l',col='red',lty=2)
lines(input_seq,m_sgasp_mle_predict@mean,type='l',col='blue',lty=2)

lines(input_seq,pred_l2$mean,col='brown',lty=2)
lines(input_seq,pred_ls_mean,col='orange',lty=2)

lines(input_all,output_all,type='p',col='black',pch=20)


legend(34,17,legend=c('N-D, PS','S-GaSP, PS','GaSP, PS',expression(L[2]),
                            'N-D, MLE','S-GaSP, MLE','GaSP, MLE','LS'),
       lty=c(1,1,1,2,2,2,2,2),col=c('green','blue','red','brown','green','blue','red','orange'),ncol=2, x.intersp=1,cex=.8)

dev.off()



pdf('post_sample_1_eg_4_box.pdf',height=4.5,width=6)

index=1
hist(m_no_discrepancy@post_sample[,index], breaks=50,
     col=rgb(0,1,0, 0.5),xlim=c(0.5,1.5),ylim=c(0,1600),xlab=expression(theta[1]),ylab='counts',main='')
hist(m_sgasp@post_sample[,index],col=rgb(0,0,1, 0.5), breaks=50,add=T)
hist(m_gasp@post_sample[,index],col=rgb(1,0,0, 0.5),breaks=50, add=T)
lines(x=m_no_discrepancy_mle@param_est[index],y=0,col=rgb(0,1,0, 0.5),cex=1.5,type='p',pch=15)
lines(x=m_gasp_mle@param_est[index],y=0,col=rgb(1,0,0, 0.5),cex=1.5,type='p',pch=17)

lines(x=m_sgasp_mle@param_est[index],y=0,col=rgb(0,0,1, 0.5),cex=1.5,type='p',pch=19)

lines(x=m_l2_record$par[index],y=0,col='brown',cex=1.5,type='p',pch=8)
lines(x=m_ls_record$par[index],y=0,col='orange',cex=1.5,type='p',pch=18)


dev.off()



pdf('post_sample_2_eg_4_box.pdf',height=4.5,width=6)
index=2
hist(m_no_discrepancy@post_sample[,index], breaks=50,
     col=rgb(0,1,0, 0.5),xlim=c(0.5,1.5),ylim=c(0,1200),xlab=expression(theta[2]),ylab='counts',main='')
hist(m_sgasp@post_sample[,index],col=rgb(0,0,1, 0.5), breaks=50,add=T)
hist(m_gasp@post_sample[,index],col=rgb(1,0,0, 0.5),breaks=50, add=T)
lines(x=m_no_discrepancy_mle@param_est[index],y=0,col=rgb(0,1,0, 0.5),cex=1.5,type='p',pch=15)
lines(x=m_gasp_mle@param_est[index],y=0,col=rgb(1,0,0, 0.5),cex=1.5,type='p',pch=17)
lines(x=m_sgasp_mle@param_est[index],y=0,col=rgb(0,0,1, 0.5),cex=1.5,type='p',pch=19)

lines(x=m_l2_record$par[index],y=0,col='brown',cex=1.5,type='p',pch=8)
lines(x=m_ls_record$par[index],y=0,col='orange',cex=1.5,type='p',pch=18)

legend('topright', legend=c('N-D, PS','S-GaSP, PS','GaSP, PS',expression(L[2]),
                           'N-D, MLE','S-GaSP, MLE','GaSP, MLE','LS'),lwd=c(4,4,4,NA,NA,NA,NA,NA),
       col=c(rgb(0,1,0, 0.5),rgb(0,0,1, 0.5),rgb(1,0,0, 0.5),'brown', rgb(0,1,0, 0.5),rgb(0,0,1, 0.5),rgb(1,0,0, 0.5),'orange')
       ,pch=c(NA,NA,NA,8,15,17,19,18),ncol=2,
       x.intersp=1,cex=.8 )
# legend('topright',legend=c("No-discrepancy","S-GaSP","GaSP"),lwd=4,
#        col=c(rgb(0,1,0, 0.5),rgb(0,0,1, 0.5),rgb(1,0,0, 0.5)))

dev.off()


# ####leave one out prediction

M=6
RMSE_all_record=matrix(0,M,8)
RMSE_fm_record=matrix(0,M,8)
avg_interval_record=matrix(0,M,8)
avg_prop_record=matrix(0,M,8)

##
avg_interval_record=matrix(0,M,8)
avg_prop_record=matrix(0,M,8)

set.seed(1)

for(i_M in 1:M){
  print(i_M)
  testing_index=(i_M-1)*2+c(1,2)



   output=as.matrix(output_mat_all[-i_M,])
   input=as.matrix(input_mat_all[-i_M])

   testing_output=t(as.matrix(output_mat_all[i_M,]))
   testing_input=as.matrix(input_all[i_M])


  theta_range=matrix(c(0.5,0.5,1.5,1.5),2,2)


  m_no_discrepancy=rcalibration(input,output,simul_type=1,math_model=Box_model_solved,
                                sd_proposal=c(0.25,0.25,1,1),
                       theta_range=theta_range,S=12000,S_0=2000,discrepancy_type = 'no-discrepancy')
  m_no_discrepancy_pred=predict(m_no_discrepancy,input_all,math_model=Box_model_solved,n_thinning=5,
                                interval_est=c(0.025, 0.975),interval_data=T)

  RMSE_fm_record[i_M,5]=sqrt(mean( (m_no_discrepancy_pred@math_model_mean[testing_index]- testing_output)^2))

  avg_interval_record[i_M,5]=mean(m_no_discrepancy_pred@interval[testing_index,2]-m_no_discrepancy_pred@interval[testing_index,1])

  avg_prop_record[i_M,5]=length(which( (testing_output> m_no_discrepancy_pred@interval[testing_index,1]) & (testing_output<m_no_discrepancy_pred@interval[testing_index,2]) ))/2



  ###gasp


  m_gasp=rcalibration(input,output,simul_type=1,math_model=Box_model_solved,
                      sd_proposal=c(0.25,0.25,1,1),
                      theta_range=theta_range,S=12000,S_0=2000,discrepancy_type = 'GaSP')


  m_gasp_pred=predict(m_gasp,input_all,math_model=Box_model_solved,n_thinning=5,
                      interval_est=c(0.025, 0.975),interval_data=T)


  RMSE_fm_record[i_M,1]=sqrt(mean( (m_gasp_pred@math_model_mean[testing_index]- testing_output)^2))

  RMSE_all_record[i_M,1]=sqrt(mean( (m_gasp_pred@mean[testing_index]- testing_output)^2))
  avg_interval_record[i_M,1]=mean(m_gasp_pred@interval[testing_index,2]-m_gasp_pred@interval[testing_index,1])

  avg_prop_record[i_M,1]=length(which( (testing_output> m_gasp_pred@interval[testing_index,1]) & (testing_output<m_gasp_pred@interval[testing_index,2]) ))/2

  ###s-gasp
  m_sgasp=rcalibration(input,output,simul_type=1,math_model=Box_model_solved,
                       sd_proposal=c(0.25,0.25,1,1),
                       theta_range=theta_range,S=12000,S_0=2000)


  model_sgasp_pred=predict(m_sgasp,input_all,math_model=Box_model_solved,n_thinning=5,
                          interval_est=c(0.025, 0.975),interval_data=T)


  RMSE_fm_record[i_M,2]=sqrt(mean( (model_sgasp_pred@math_model_mean[testing_index]- testing_output)^2))

  RMSE_all_record[i_M,2]=sqrt(mean( (model_sgasp_pred@mean[testing_index]- testing_output)^2))

  avg_interval_record[i_M,2]=mean(model_sgasp_pred@interval[testing_index,2]-model_sgasp_pred@interval[testing_index,1])

  avg_prop_record[i_M,2]=length(which( (testing_output> model_sgasp_pred@interval[testing_index,1]) & (testing_output<model_sgasp_pred@interval[testing_index,2]) ))/2



  ###two steps
  output_vec=as.vector(t(output))
  input_vec=as.vector(t(matrix(input,length(input),2)))
  ##l2
  m_l2=rgasp(input_vec,output_vec,nugget.est=T)
  pred_l2=predict(m_l2,input_all,interval_data=T)


  avg_interval_record[i_M,3]=mean(pred_l2$upper95[testing_index]-pred_l2$lower95[testing_index])

  avg_prop_record[i_M,3]=length(which( (testing_output> pred_l2$lower95[testing_index]) & (testing_output<pred_l2$upper95[testing_index]) ))/2

  input_seq=as.matrix(seq(input_all[1],input_all[12],.1))
  pred_l2_input_seq=predict(m_l2,input_seq)


  m_l2_record=try(optim(c(1,1), L2_calibration, testing_input_here=input_seq, pred_testing_output=pred_l2_input_seq$mean,
                        method="L-BFGS-B",lower=theta_range[,1],upper=theta_range[,2]),silent = T)

  pred_l2_cm=Box_model_solved(input_all,m_l2_record$par)

  RMSE_fm_record[i_M,3]=sqrt(mean( (pred_l2_cm[testing_index]-testing_output)^2))
  RMSE_all_record[i_M,3]=sqrt(mean( (pred_l2$mean[testing_index]-testing_output)^2))


  ##ls
  m_ls_record=try(optim(c(1,1), LS_calibration,input_here=input_vec,output_here=output_vec,
                        method="L-BFGS-B",lower=theta_range[,1],upper=theta_range[,2]),silent = T)

  #record_par_ls[i_M,]=m_ls_record$par

  output_tilde_ls=output_vec-Box_model_solved(input_vec,m_ls_record$par)

  m_ls=rgasp(input_vec,output_tilde_ls,nugget.est=T)
  pred_ls_cm_mean=Box_model_solved(input_all,m_ls_record$par)
  m_ls_pred=predict(m_ls,input_all)
  pred_ls_mean=m_ls_pred$mean+pred_ls_cm_mean

  RMSE_fm_record[i_M,4]=sqrt(mean( (pred_ls_cm_mean[testing_index]-testing_output)^2))
  RMSE_all_record[i_M,4]=sqrt(mean( (pred_ls_mean[testing_index]-testing_output)^2))


  pred_ls_lower95=m_ls_pred$lower95+pred_ls_cm_mean
  pred_ls_upper95=m_ls_pred$upper95+pred_ls_cm_mean

  avg_interval_record[i_M,4]=mean(pred_ls_upper95[testing_index]-pred_ls_lower95[testing_index])

  avg_prop_record[i_M,4]=length(which( (testing_output> pred_ls_lower95[testing_index]) & (testing_output<pred_ls_upper95[testing_index]) ))/2

  ###2021 add optimize
  #unloadNamespace('RobustCalibration')
  #library(RobustCalibration)


  m_gasp_mle=rcalibration(design=input, observations=output,p_theta=2,simul_type=1,math_model=Box_model_solved,
                   theta_range=theta_range,discrepancy_type = 'GaSP',method='mle',num_initial_values=5)

  m_sgasp_mle=rcalibration(design=input, observations=output,p_theta=2,simul_type=1,math_model=Box_model_solved,
                          theta_range=theta_range,discrepancy_type = 'S-GaSP',method='mle',num_initial_values=5)
  m_no_discrepancy_mle=rcalibration(design=input, observations=output,simul_type=1,math_model=Box_model_solved,
                                    theta_range=theta_range,discrepancy_type = 'no-discrepancy',method='mle',num_initial_values=5)



  m_gasp_mle_predict=predict(m_gasp_mle,input_all,math_model=Box_model_solved,
                                       interval_est=c(0.025, 0.975),interval_data=T)  ##no mean, maybe set mean to be math model mean?

  m_sgasp_mle_predict=predict(m_sgasp_mle,input_all,math_model=Box_model_solved,
                             interval_est=c(0.025, 0.975),interval_data=T)  ##no mean, maybe set mean to be math model mean?

  m_no_discrepancy_mle_predict=predict(m_no_discrepancy_mle,input_all,math_model=Box_model_solved,
                                       interval_est=c(0.025, 0.975),interval_data=T)  ##no mean, maybe set mean to be math model mean?


  ###gasp, mle
  RMSE_fm_record[i_M,6]=sqrt(mean( (m_gasp_mle_predict@math_model_mean[testing_index]- testing_output)^2))

  RMSE_all_record[i_M,6]=sqrt(mean( (m_gasp_mle_predict@mean[testing_index]- testing_output)^2))

  avg_interval_record[i_M,6]=mean(m_gasp_mle_predict@interval[testing_index,2]-m_gasp_mle_predict@interval[testing_index,1])

  avg_prop_record[i_M,6]=length(which( (testing_output> m_gasp_mle_predict@interval[testing_index,1]) & (testing_output<m_gasp_mle_predict@interval[testing_index,2]) ))/2


  ##sgasp, mle
  RMSE_fm_record[i_M,7]=sqrt(mean( (m_sgasp_mle_predict@math_model_mean[testing_index]- testing_output)^2))

  RMSE_all_record[i_M,7]=sqrt(mean( (m_sgasp_mle_predict@mean[testing_index]- testing_output)^2))

  avg_interval_record[i_M,7]=mean(m_sgasp_mle_predict@interval[testing_index,2]-m_sgasp_mle_predict@interval[testing_index,1])

  avg_prop_record[i_M,7]=length(which( (testing_output> m_sgasp_mle_predict@interval[testing_index,1]) & (testing_output<m_sgasp_mle_predict@interval[testing_index,2]) ))/2


  ##sgasp, no-discrepancy
  RMSE_fm_record[i_M,8]=sqrt(mean( (m_no_discrepancy_mle_predict@math_model_mean[testing_index]- testing_output)^2))

  RMSE_all_record[i_M,8]=sqrt(mean( (m_no_discrepancy_mle_predict@mean[testing_index]- testing_output)^2))

  avg_interval_record[i_M,8]=mean(m_no_discrepancy_mle_predict@interval[testing_index,2]-m_no_discrepancy_mle_predict@interval[testing_index,1])

  avg_prop_record[i_M,8]=length(which( (testing_output> m_no_discrepancy_mle_predict@interval[testing_index,1]) & (testing_output<m_no_discrepancy_mle_predict@interval[testing_index,2]) ))/2

  #plot(m_no_discrepancy@post_sample


  print(RMSE_fm_record[i_M,])
  print(RMSE_all_record[i_M,])

}
# 
# 
colMeans(RMSE_fm_record)
colMeans(RMSE_all_record)
colMeans(avg_interval_record)
colMeans(avg_prop_record)


