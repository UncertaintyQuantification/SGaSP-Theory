###Ion Channel model
library(expm)
library(RobustCalibration)
library(RobustGaSP)

input_real=as.matrix(c(-1.7148,-1.0217,-0.5798,-0.1985, 0.0953, 0.3507,0.5988, 0.8242,1.0367,1.2355,1.4303,
          1.6214,1.8050,1.9851,2.1633,2.3380,2.5096,2.6797,2.8495))
output_real= as.matrix(c(0.0538,0.0878,0.1211,0.1356,0.1184,0.0900,0.0649,0.0471,0.0338,0.0244,
0.0175,0.0118,0.0089,0.0060,0.0045,0.0027,0.0024,0.0018,0.0013))

plot(input_real,output_real)



##x is the logarithm of time
ion_channel_model<-function(x,theta){
  A=matrix(0,4,4)
  A[1,1]=-theta[2]-theta[3]
  A[1,2]=theta[1]
  A[2,1]=theta[2]
  A[2,2]=-theta[1]-theta[2]
  A[2,3]=theta[1]
  A[3,2]=theta[2]
  A[3,3]=-theta[1]-theta[2]
  A[3,4]=theta[1]
  A[4,3]=theta[2]
  A[4,4]=-theta[1]
  n_x=length(x)
  res=rep(0,n_x)
  for(i_x in 1:n_x){
    res[i_x]=t(c(1,0,0,0))%*%expm(exp(x[i_x])*A)%*%c(0,0,0,1)
  }
  res
}


LS_calibration<-function(theta,input_here,output_here){
  mean((ion_channel_model(input_here,theta)-output_here)^2)
}

L2_calibration<-function(theta,testing_input_here,pred_testing_output){
  mean((ion_channel_model(testing_input_here,theta)-pred_testing_output)^2)
}


theta_range=matrix(c(0,0,0,10,10,10),3,2)


M=length(output_real)
RMSE_all_record=matrix(0,M,8) 
RMSE_fm_record=matrix(0,M,8) 

abs_err_all_record=matrix(0,M,8) 
abs_err_fm_record=matrix(0,M,8) 


##
avg_interval_record=matrix(0,M,8) 
avg_prop_record=matrix(0,M,8) 


record_time=system.time(
  for(i_M in 1:M){
    print(i_M)
    set.seed(i_M)
    testing_index=i_M
    
    
    #
    m_no_discrepancy=rcalibration(design=input_real[-testing_index], observations=output_real[-testing_index],simul_type=1,
                                  theta_range = theta_range,math_model = ion_channel_model,
                                  sd_proposal=c(rep(0.03,3),0.3,0.3),
                                  discrepancy_type='no-discrepancy',S=12000,S_0=2000)
    m_no_discrepancy_pred=predict(m_no_discrepancy,input_real,interval_est=c(0.025,0.975),
                                  math_model=ion_channel_model,n_thinning=10,interval_data=T)
    
    
    RMSE_fm_record[i_M,5]=sqrt(mean( (m_no_discrepancy_pred@math_model_mean[testing_index]- output_real[testing_index])^2))
    
    
    avg_interval_record[i_M,5]=mean(m_no_discrepancy_pred@interval[testing_index,2]-m_no_discrepancy_pred@interval[testing_index,1])
    
    avg_prop_record[i_M,5]=length(which( (output_real[testing_index]> m_no_discrepancy_pred@interval[testing_index,1]) & (output_real[testing_index]<m_no_discrepancy_pred@interval[testing_index,2]) ))/1
    
    ##GaSP
    m_GaSP=rcalibration(design=input_real[-testing_index], observations=output_real[-testing_index],p_theta=3,simul_type=1,
                        theta_range = theta_range,math_model = ion_channel_model,
                        sd_proposal=c(rep(0.03,3),0.3,0.3),
                        discrepancy_type='GaSP',S=12000,S_0=2000)
    
    m_GaSP_pred=predict(m_GaSP,input_real,interval_est=c(0.025,0.975),
                        math_model=ion_channel_model,n_thinning=10,interval_data=T)
    
    
    RMSE_fm_record[i_M,1]=sqrt(mean( (m_GaSP_pred@math_model_mean[testing_index]- output_real[testing_index])^2))
    
    RMSE_all_record[i_M,1]=sqrt(mean( (m_GaSP_pred@mean[testing_index]- output_real[testing_index])^2))
    avg_interval_record[i_M,1]=mean(m_GaSP_pred@interval[testing_index,2]-m_GaSP_pred@interval[testing_index,1])
    
    avg_prop_record[i_M,1]=length(which( (output_real[testing_index]> m_GaSP_pred@interval[testing_index,1]) & (output_real[testing_index]<m_GaSP_pred@interval[testing_index,2]) ))/1
    
    
    ##SGaSP
    m_SGaSP=rcalibration(design=input_real[-testing_index], observations=output_real[-testing_index],p_theta=3,simul_type=1,
                         theta_range = theta_range,math_model = ion_channel_model,
                         sd_proposal=c(rep(0.03,3),0.3,0.3),
                         discrepancy_type='S-GaSP',S=12000,S_0=2000)

    

    

    m_SGaSP_pred=predict(m_SGaSP,input_real,interval_est=c(0.025,0.975),
                         math_model=ion_channel_model,n_thinning=10,interval_data=T)
    


    RMSE_fm_record[i_M,2]=sqrt(mean( (m_SGaSP_pred@math_model_mean[testing_index]- output_real[testing_index])^2))
    
    RMSE_all_record[i_M,2]=sqrt(mean( (m_SGaSP_pred@mean[testing_index]- output_real[testing_index])^2))
    
    avg_interval_record[i_M,2]=mean(m_SGaSP_pred@interval[testing_index,2]-m_SGaSP_pred@interval[testing_index,1])
    
    avg_prop_record[i_M,2]=length(which( (output_real[testing_index]> m_SGaSP_pred@interval[testing_index,1]) & (output_real[testing_index]<m_SGaSP_pred@interval[testing_index,2]) ))/1
    
    
    
    ###two steps
    ##l2

    m_l2=rgasp(input_real[-testing_index],output_real[-testing_index],nugget.est=T)
    pred_l2=predict(m_l2,input_real,interval_data=T)
    
    
    avg_interval_record[i_M,3]=mean(pred_l2$upper95[testing_index]-pred_l2$lower95[testing_index])
    
    avg_prop_record[i_M,3]=length(which( (output_real[testing_index]> pred_l2$lower95[testing_index]) & (output_real[testing_index]<pred_l2$upper95[testing_index]) ))/1
    
    input_seq=as.matrix(seq(-2,3,5/499))
    pred_l2_input_seq=predict(m_l2,input_seq)
    
    
    m_l2_record=try(optim(rowMeans(theta_range), L2_calibration, testing_input_here=input_seq, pred_testing_output=pred_l2_input_seq$mean, 
                          method="L-BFGS-B",lower=theta_range[,1],upper=theta_range[,2]),silent = T)
    
    pred_l2_cm=ion_channel_model(input_real,m_l2_record$par)
    
    RMSE_fm_record[i_M,3]=sqrt(mean( (pred_l2_cm[testing_index]-output_real[testing_index])^2))
    RMSE_all_record[i_M,3]=sqrt(mean( (pred_l2$mean[testing_index]-output_real[testing_index])^2))
    
    
    ##ls
    m_ls_record=try(optim(rowMeans(theta_range), LS_calibration,input_here=input_real[-testing_index],output_here=output_real[-testing_index],
                          method="L-BFGS-B",lower=theta_range[,1],upper=theta_range[,2]),silent = T)
    
    #record_par_ls[i_M,]=m_ls_record$par
    
    output_tilde_ls=output_real[-testing_index]-ion_channel_model(input_real[-testing_index],m_ls_record$par)
    
    m_ls=rgasp(input_real[-testing_index],output_tilde_ls,nugget.est=T)
    
    pred_ls_cm_mean=ion_channel_model(input_real,m_ls_record$par)
    m_ls_pred=predict(m_ls,input_real)
    pred_ls_mean=m_ls_pred$mean+pred_ls_cm_mean
    
    RMSE_fm_record[i_M,4]=sqrt(mean( (pred_ls_cm_mean[testing_index]-output_real[testing_index])^2))
    RMSE_all_record[i_M,4]=sqrt(mean( (pred_ls_mean[testing_index]-output_real[testing_index])^2))
    
    
    pred_ls_lower95=m_ls_pred$lower95+pred_ls_cm_mean
    pred_ls_upper95=m_ls_pred$upper95+pred_ls_cm_mean
    
    avg_interval_record[i_M,4]=mean(pred_ls_upper95[testing_index]-pred_ls_lower95[testing_index])
    
    avg_prop_record[i_M,4]=length(which( (output_real[testing_index]> pred_ls_lower95[testing_index]) & (output_real[testing_index]<pred_ls_upper95[testing_index]) ))/1
    

    ##added MLE approaches
    
    
    ###2021 add optimize
    ##
    #unloadNamespace('RobustCalibration')
    #library(RobustCalibration)
    
    

    m_gasp_mle=rcalibration(design=input_real[-testing_index], observations=output_real[-testing_index],simul_type=1,math_model=ion_channel_model,
                            theta_range=theta_range,discrepancy_type = 'GaSP',method='mle',num_initial_values=5)
    
    m_sgasp_mle=rcalibration(design=input_real[-testing_index], observations=output_real[-testing_index],simul_type=1,math_model=ion_channel_model,
                             theta_range=theta_range,discrepancy_type = 'S-GaSP',method='mle',num_initial_values=5)
    m_no_discrepancy_mle=rcalibration(design=input_real[-testing_index], observations=output_real[-testing_index],simul_type=1,math_model=ion_channel_model,
                                      theta_range=theta_range,discrepancy_type = 'no-discrepancy',method='mle',num_initial_values=5)
    
    
    
    m_gasp_mle_predict=predict(m_gasp_mle,input_real,math_model=ion_channel_model,
                               interval_est=c(0.025, 0.975),interval_data=T)  ##no mean, maybe set mean to be math model mean?
    
    m_sgasp_mle_predict=predict(m_sgasp_mle,input_real,math_model=ion_channel_model,
                                interval_est=c(0.025, 0.975),interval_data=T)  ##no mean, maybe set mean to be math model mean?
    
    m_no_discrepancy_mle_predict=predict(m_no_discrepancy_mle,input_real,math_model=ion_channel_model,
                                         interval_est=c(0.025, 0.975),interval_data=T)  ##no mean, maybe set mean to be math model mean?
    
    
    ###gasp, mle
    RMSE_fm_record[i_M,6]=sqrt(mean( (m_gasp_mle_predict@math_model_mean[testing_index]-  output_real[testing_index])^2))
    
    RMSE_all_record[i_M,6]=sqrt(mean( (m_gasp_mle_predict@mean[testing_index]-  output_real[testing_index])^2))
    
    avg_interval_record[i_M,6]=mean(m_gasp_mle_predict@interval[testing_index,2]-m_gasp_mle_predict@interval[testing_index,1])
    
    avg_prop_record[i_M,6]=length(which( ( output_real[testing_index]> m_gasp_mle_predict@interval[testing_index,1]) & ( output_real[testing_index]<m_gasp_mle_predict@interval[testing_index,2]) ))
    
    
    

    
    
    ##sgasp, mle
    RMSE_fm_record[i_M,7]=sqrt(mean( (m_sgasp_mle_predict@math_model_mean[testing_index]- output_real[testing_index])^2))
    
    RMSE_all_record[i_M,7]=sqrt(mean( (m_sgasp_mle_predict@mean[testing_index]- output_real[testing_index])^2))
    
    avg_interval_record[i_M,7]=mean(m_sgasp_mle_predict@interval[testing_index,2]-m_sgasp_mle_predict@interval[testing_index,1])
    
    avg_prop_record[i_M,7]=length(which( (output_real[testing_index]> m_sgasp_mle_predict@interval[testing_index,1]) & (output_real[testing_index]<m_sgasp_mle_predict@interval[testing_index,2]) ))/1
    
    
    ##sgasp, no-discrepancy
    RMSE_fm_record[i_M,8]=sqrt(mean( (m_no_discrepancy_mle_predict@math_model_mean[testing_index]- output_real[testing_index])^2))
    
    RMSE_all_record[i_M,8]=sqrt(mean( (m_no_discrepancy_mle_predict@mean[testing_index]- output_real[testing_index])^2))
    
    avg_interval_record[i_M,8]=mean(m_no_discrepancy_mle_predict@interval[testing_index,2]-m_no_discrepancy_mle_predict@interval[testing_index,1])
    
    avg_prop_record[i_M,8]=length(which( (output_real[testing_index]> m_no_discrepancy_mle_predict@interval[testing_index,1]) & (output_real[testing_index]<m_no_discrepancy_mle_predict@interval[testing_index,2]) ))/1
    
    #plot(m_no_discrepancy@post_sample
    
    
    ##  
    print(RMSE_fm_record[i_M,])
    print(RMSE_all_record[i_M,])
    print(avg_prop_record[i_M,])
    print(avg_interval_record[i_M,])
    
  }
)




colMeans(RMSE_fm_record)
colMeans(RMSE_all_record)
colMeans(avg_interval_record)
colMeans(avg_prop_record)

# GaSP PS, S-GaSP PS, L2, LS, ND PS, GaSP MLE, S-GaSP MLE,  ND PS
# > colMeans(RMSE_fm_record)
# [1] 0.011297883 0.005624674 0.004836469 0.006494430 0.006452842 0.017456522 0.005763965 0.006508024
# > colMeans(RMSE_all_record)
# [1] 0.0017655382 0.0015248960 0.0008260662 0.0021654446 0.0000000000 0.0019731873 0.0013597069
# [8] 0.0065080241
# > colMeans(avg_interval_record)
# [1] 0.005646011 0.009580026 0.009446482 0.005797439 0.025802990 0.003564241 0.007101974 0.022364874
# > colMeans(avg_prop_record)
# [1] 0.9473684 1.0000000 0.8421053 0.7894737 0.8421053 0.8421053 0.7368421 0.8421053


#
m_no_discrepancy=rcalibration(design=input_real, observations=output_real,simul_type=1,
                              theta_range = theta_range,math_model = ion_channel_model,sd_proposal=c(rep(0.03,3),0.3,0.3),
                              discrepancy_type='no-discrepancy',S=12000,S_0=2000)

m_GaSP=rcalibration(design=input_real, observations=output_real,simul_type=1,
                              theta_range = theta_range,math_model = ion_channel_model,sd_proposal=c(rep(0.03,3),0.3,0.3),
                              discrepancy_type='GaSP',S=12000,S_0=2000)

m_SGaSP=rcalibration(design=input_real, observations=output_real,simul_type=1,
                    theta_range = theta_range,math_model = ion_channel_model,sd_proposal=c(rep(0.03,3),0.3,0.3),
                    discrepancy_type='S-GaSP',S=12000,S_0=2000)


input_seq=as.matrix(seq(-2,3,5/499))


m_no_discrepancy_pred=predict(m_no_discrepancy,input_seq,interval_est=c(0.025,0.975),
                              math_model=ion_channel_model,n_thinning=10)

m_GaSP_pred=predict(m_GaSP,input_seq,interval_est=c(0.025,0.975),
                    math_model=ion_channel_model,n_thinning=10)
m_SGaSP_pred=predict(m_SGaSP,input_seq,interval_est=c(0.025,0.975),
                     math_model=ion_channel_model,n_thinning=10)


###MLE
m_no_discrepancy_mle=rcalibration(design=input_real, observations=output_real,simul_type=1,
                              theta_range = theta_range,math_model = ion_channel_model,sd_proposal=c(rep(0.03,3),0.3,0.3),
                              discrepancy_type='no-discrepancy',method='mle',num_initial_values = 5)

m_GaSP_mle=rcalibration(design=input_real, observations=output_real,simul_type=1,
                    theta_range = theta_range,math_model = ion_channel_model,sd_proposal=c(rep(0.03,3),0.3,0.3),
                    discrepancy_type='GaSP',method='mle',num_initial_values = 5)

m_SGaSP_mle=rcalibration(design=input_real, observations=output_real,simul_type=1,
                     theta_range = theta_range,math_model = ion_channel_model,sd_proposal=c(rep(0.03,3),0.3,0.3),
                     discrepancy_type='S-GaSP',method='mle',num_initial_values = 5)


m_gasp_mle_predict=predict(m_gasp_mle,input_seq,math_model=ion_channel_model,
                           interval_est=c(0.025, 0.975),interval_data=T)  

m_sgasp_mle_predict=predict(m_sgasp_mle,input_seq,math_model=ion_channel_model,
                            interval_est=c(0.025, 0.975),interval_data=T)  

m_no_discrepancy_mle_predict=predict(m_no_discrepancy_mle,input_seq,math_model=ion_channel_model,
                                     interval_est=c(0.025, 0.975),interval_data=T)  



###two steps
##l2

m_l2=rgasp(input_real,output_real,nugget.est=T)
pred_l2_input_seq=predict(m_l2,input_seq,interval_data=T)

m_l2_record=try(optim(rowMeans(theta_range), L2_calibration, testing_input_here=input_seq, pred_testing_output=pred_l2_input_seq$mean, 
                      method="L-BFGS-B",lower=theta_range[,1],upper=theta_range[,2]),silent = T)

pred_l2_cm=ion_channel_model(input_seq,m_l2_record$par)

##ls
m_ls_record=try(optim(rowMeans(theta_range), LS_calibration,input_here=input_real,output_here=output_real,
                      method="L-BFGS-B",lower=theta_range[,1],upper=theta_range[,2]),silent = T)


output_tilde_ls=output_real-ion_channel_model(input_real,m_ls_record$par)

m_ls=rgasp(input_real,output_tilde_ls,nugget.est=T)

pred_ls_cm_mean=ion_channel_model(input_seq,m_ls_record$par)
m_ls_pred=predict(m_ls,input_seq)
pred_ls_mean=m_ls_pred$mean+pred_ls_cm_mean



pred_ls_lower95=m_ls_pred$lower95+pred_ls_cm_mean
pred_ls_upper95=m_ls_pred$upper95+pred_ls_cm_mean





pdf('plot_calibration_eg_4_ion_channel.pdf',height=4.33,width=6.5)

plot(input_real,output_real,type='p',pch=20,xlab='x',ylab='y',
     main='Calibrated computer model',mgp=c(2.5,1,0))
polygon( c(input_seq,rev(input_seq)),c(m_GaSP_pred@interval[,1], rev(m_GaSP_pred@interval[,2])),col =  "grey80", border = FALSE)
lines(input_real,output_real,type='p',pch=20)
lines(input_seq,m_no_discrepancy_pred@math_model_mean,col='green',lty=1)
lines(input_seq,m_GaSP_pred@math_model_mean,col='red',lty=1)
lines(input_seq,m_SGaSP_pred@math_model_mean,col='blue',lty=1)
lines(input_seq,m_no_discrepancy_mle_predict@math_model_mean,col='green',lty=2)
lines(input_seq,m_gasp_mle_predict@math_model_mean,col='red',lty=2)
lines(input_seq,m_sgasp_mle_predict@math_model_mean,col='blue',lty=2)
lines(input_seq,pred_l2_cm,col='brown',lty=2)
lines(input_seq,pred_ls_cm_mean,col='orange',lty=2)

legend("topright", legend=c('N-D, PS','S-GaSP, PS','GaSP, PS',expression(L[2]),
                          'N-D, MLE','S-GaSP, MLE','GaSP, MLE','LS'),
       lty=c(1,1,1,2,2,2,2,2),col=c('green','blue','red','brown','green','blue','red','orange'),ncol=2, x.intersp=1,cex=.8)

dev.off()

pdf('plot_pred_eg_4_ion_channel.pdf',height=4.33,width=6.5)

plot(input_real,output_real,type='p',pch=20,xlab='x',ylab='y',
     main='Calibrated computer model+discrepancy',mgp=c(2.5,1,0))
polygon( c(input_seq,rev(input_seq)),c(m_SGaSP_pred@interval[,1], rev(m_SGaSP_pred@interval[,2])),col =  "grey80", 
         border = FALSE)
lines(input_real,output_real,type='p',pch=20)
lines(input_seq,m_GaSP_pred@mean,col='red',lty=1)
lines(input_seq,m_SGaSP_pred@mean,col='blue',lty=1)
lines(input_seq,m_gasp_mle_predict@mean,col='red',lty=2)
lines(input_seq,m_sgasp_mle_predict@mean,col='blue',lty=2)
lines(input_seq,pred_l2_input_seq$mean,col='brown',lty=2)
lines(input_seq,pred_ls_mean,col='orange',lty=2)

legend("topright", legend=c('N-D, PS','S-GaSP, PS','GaSP, PS',expression(L[2]),
                            'N-D, MLE','S-GaSP, MLE','GaSP, MLE','LS'),
       lty=c(1,1,1,2,2,2,2,2),col=c('green','blue','red','brown','green','blue','red','orange'),ncol=2, x.intersp=1,cex=.8)

dev.off()





pdf('post_sample_1_2_eg_5_ion_channel.pdf',height=5,width=6.3)
plot(m_no_discrepancy@post_sample[seq(1,1000)*10,1],m_no_discrepancy@post_sample[seq(1,1000)*10,2],type='p',
     pch=15,col='green',xlim=c(theta_range[1,1],theta_range[1,2]),ylim=c(theta_range[2,1],theta_range[2,2]),
     cex=0.6,xlab=expression(theta[1]),ylab=expression(theta[2]),mgp=c(2.5,1,0) )
lines(m_SGaSP@post_sample[seq(1,1000)*10,1],m_SGaSP@post_sample[seq(1,1000)*10,2],type='p',
      pch=20,col='blue',cex=0.6)
lines(m_GaSP@post_sample[seq(1,1000)*10,1],m_GaSP@post_sample[seq(1,1000)*10,2],type='p',
      pch=17,col='red',cex=0.6)
# legend("topleft",  legend=c('N-D, PS','S-GaSP, PS','GaSP, PS',
#                             'N-D, MLE','S-GaSP, MLE','GaSP, MLE'), pch=c(15,20,17,NA,NA,NA),
#        lty=c(NA,NA,NA,1,2,3),col=c('green','blue','red','green','blue','red'),ncol=2,
#        x.intersp=1,cex=.8)
abline(v=m_no_discrepancy_mle@param_est[1],h=m_no_discrepancy_mle@param_est[2],lty=1,col='green')
abline(v=m_GaSP_mle@param_est[1],h=m_GaSP_mle@param_est[2],lty=2,col='red')
abline(v=m_SGaSP_mle@param_est[1],h=m_SGaSP_mle@param_est[2],lty=3,col='blue')
abline(v=m_l2_record$par[1],h=m_l2_record$par[2],lty=4,col='brown')
abline(v=m_ls_record$par[1],h=m_ls_record$par[2],lty=5,col='orange')

dev.off()
pdf('post_sample_1_3_eg_5_ion_channel.pdf',height=5,width=6.3)
plot(m_no_discrepancy@post_sample[seq(1,1000)*10,1],m_no_discrepancy@post_sample[seq(1,1000)*10,3],type='p',
     pch=15,col='green',xlim=c(theta_range[1,1],theta_range[1,2]),ylim=c(theta_range[3,1],theta_range[3,2]),cex=0.6,
     xlab=expression(theta[1]),ylab=expression(theta[3]),mgp=c(2.5,1,0) )
lines(m_SGaSP@post_sample[seq(1,1000)*10,1],m_SGaSP@post_sample[seq(1,1000)*10,3],type='p',
      pch=20,col='blue',cex=0.6)
lines(m_GaSP@post_sample[seq(1,1000)*10,1],m_GaSP@post_sample[seq(1,1000)*10,3],type='p',
      pch=17,col='red',cex=0.6)
abline(v=m_no_discrepancy_mle@param_est[1],h=m_no_discrepancy_mle@param_est[3],lty=1,col='green')
abline(v=m_GaSP_mle@param_est[1],h=m_GaSP_mle@param_est[3],lty=2,col='red')
abline(v=m_SGaSP_mle@param_est[1],h=m_SGaSP_mle@param_est[3],lty=3,col='blue')
abline(v=m_l2_record$par[1],h=m_l2_record$par[3],lty=4,col='brown')
abline(v=m_ls_record$par[1],h=m_ls_record$par[3],lty=5,col='orange')


dev.off()


pdf('post_sample_2_3_eg_5_ion_channel.pdf',height=5,width=6.3)

plot(m_no_discrepancy@post_sample[seq(1,1000)*10,2],m_no_discrepancy@post_sample[seq(1,1000)*10,3],type='p',
     pch=15,col='green',xlim=c(theta_range[2,1],theta_range[2,2]),ylim=c(theta_range[3,1],theta_range[3,2]),cex=0.6,
     xlab=expression(theta[2]),ylab=expression(theta[3]),mgp=c(2.5,1,0) )
lines(m_SGaSP@post_sample[seq(1,1000)*10,2],m_SGaSP@post_sample[seq(1,1000)*10,3],type='p',
      pch=20,col='blue',cex=0.6)
lines(m_GaSP@post_sample[seq(1,1000)*10,2],m_GaSP@post_sample[seq(1,1000)*10,3],type='p',
      pch=17,col='red',cex=0.6)
abline(v=m_no_discrepancy_mle@param_est[2],h=m_no_discrepancy_mle@param_est[3],lty=1,col='green')
abline(v=m_GaSP_mle@param_est[2],h=m_GaSP_mle@param_est[3],lty=2,col='red')
abline(v=m_SGaSP_mle@param_est[2],h=m_SGaSP_mle@param_est[3],lty=3,col='blue')
abline(v=m_l2_record$par[2],h=m_l2_record$par[3],lty=4,col='brown')
abline(v=m_ls_record$par[2],h=m_ls_record$par[3],lty=5,col='orange')
legend("topright",  legend=c('N-D, PS','S-GaSP, PS','GaSP, PS',expression(L[2]),
                             'N-D, MLE','S-GaSP, MLE','GaSP, MLE','LS'), pch=c(15,20,17,NA,NA,NA,NA,NA),
       lty=c(NA,NA,NA,4,1,2,3,5),col=c('green','blue','red','brown','green','blue','red','orange'),ncol=2,
       x.intersp=1,cex=.8)

dev.off()

