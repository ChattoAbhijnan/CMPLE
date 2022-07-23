log_likelihood<- function(para){


  
  alpha<- array(0,dim=c(q,p+1))
  delta<- array(0,c(q,q,p+1))
  
  
  place1=0
  for (i in 1:q) {
    
    alpha[i,]<- para[(place1+1):(place1+p+1)]
    place1=place1+p+1
    
    
  }
  
  place2<- q*(p+1)
  for (j  in 1:(q-1) ){
    for (k in (j+1):q) {
      
      delta[j,k,]<- para[(place2+1):(place2+p+1)]
      place2=place2+p+1
    }
  }
  
  
  Sig<- array(0,dim=c(n,q))
  for (j in 1:q) {
    Sig[,j]<- exp (X %*% alpha[j,])
    
  }
  
  #Rho calculation 
  Rho<- array(0,c(n,q,q))
  
  for (j  in 1:(q-1) ){
    for (k in (j+1):q) {
      Rho[,j,k]<- 1- 2/(1+ exp(X %*% delta[j,k,]))
      
    }
    
  }
  
  for (i  in 1:n) {
    Rho[i,,]<- Rho[i,,]+ t(Rho[i,,])
    diag(Rho[i,,])<- 1
    
  }
  Rho[Rho<= -0.99]= -0.99
  Rho[Rho>= 0.999]= 0.999
  
  log_lik<- array(0,c(n,q,q))
  for (i in 1:n) {
    for (j in 1:(q-1)) {
      for (k in (j+1):q) {
        
        log_lik[i,j,k]<- (-0.5)*(2*(X[i,]%*%alpha[j,]) + 2*(X[i,]%*%alpha[k,]) + log(1- (1-2/(1+exp(X[i,]%*%delta[j,k,])))**2)+ (1/(1- (1-2/(1+exp(X[i,]%*%delta[j,k,])))**2))* ((Y[i,j]**2/(exp(2*X[i,]%*%alpha[j,]))) + (Y[i,k]**2/(exp(2*X[i,]%*%alpha[k,])))- ((2*(1-2/(1+exp(X[i,]%*%delta[j,k,])))*Y[i,j]*Y[i,k])/exp((X[i,]%*%alpha[j,])+ (X[i,]%*%alpha[k,])) )))
      }
      
    }
    
  }
  
  sum_log_lik<- apply(log_lik, c(1),sum)
  
  
  return(sum_log_lik)
  
  
}

u_theta<- function(para){

  
  alpha<- array(0,dim=c(q,p+1))
  delta<- array(0,c(q,q,p+1))
  
  
  place1=0
  for (i in 1:q) {
    
    alpha[i,]<- para[(place1+1):(place1+p+1)]
    place1=place1+p+1
    
    
  }
  
  place2<- q*(p+1)
  for (j  in 1:(q-1) ){
    for (k in (j+1):q) {
      
      delta[j,k,]<- para[(place2+1):(place2+p+1)]
      place2=place2+p+1
    }
  }
  
  
  Sig<- array(0,dim=c(n,q))
  for (j in 1:q) {
    Sig[,j]<- exp (X %*% alpha[j,])
    
  }
  
  #Rho calculation 
  Rho<- array(0,c(n,q,q))
  
  for (j  in 1:(q-1) ){
    for (k in (j+1):q) {
      Rho[,j,k]<- 1- 2/(1+ exp(X %*% delta[j,k,]))
      
    }
    
  }
  
  for (i  in 1:n) {
    Rho[i,,]<- Rho[i,,]+ t(Rho[i,,])
    diag(Rho[i,,])<- 1
    
  }
  ####
  
  
  u_alpha<- array(0,c(n,q,q,p+1))
  
  
  
  for (i in 1:n) {
    
    for (j in 1: q)
    {
      for (k in 1:q) {
        if(k==j) 
          
        {
          u_alpha[i,j,k,]<- (-(q-1))*X[i,] 
        }
        else{
          
          #print(paste("alpha[i,j,k,]",j,k))
          u_alpha[i,j,k,]<-((Y[i,j]**2/((1-Rho[i,j,k]**2)*Sig[i,j]**2)) - ((Rho[i,j,k]*Y[i,j]*Y[i,k])/((1-Rho[i,j,k]**2)*Sig[i,k]*Sig[i,j])))*X[i,]
          #(-(q-1) + (Y[i,j]**2/((Sig[i,j]**2)* (1-Rho[i,j,k]**2))) -(Rho[i,j,k]*Y[i,j]*Y[i,k]/(Sig[i,j]*Sig[i,k])))*X[i,]
        }
        
      }
      
    }
  }
  
  u_alpha_2<- apply(u_alpha,c(2,4),sum)
  
  
  u_delta<- array(0,c(n,q,q,p+1))
  
  
  for(i in 1: n){
    for (j in 1: (q-1))
    {
      for (k in (j+1):q) {
        # print(paste("del_delta[i,j,k,]",j,k))
        
        u_delta[i,j,k,]<- (0.5*Rho[i,j,k] - (0.5* ((Y[i,j]**2/Sig[i,j]**2)+ (Y[i,k]**2/Sig[i,k]**2))*(Rho[i,j,k]/(1-Rho[i,j,k]**2)))+
                             ((Y[i,j]*Y[i,k]/(Sig[i,j]*Sig[i,k]))*(0.5+ Rho[i,j,k]**2/(1-Rho[i,j,k]**2)))) *X[i,]  
        #((Y[i,j]*Y[i,k])/(Sig[i,j]*Sig[i,k])) *(0.5+ (Rho[i,j,k]**2/(1-Rho[i,j,k]**2))))*X[i,]
        
      }
    }
  }
  u_delta_2<- apply(u_delta,c(2,3,4),sum)
  
  
  
  
  ####
  
  
  all_u_alpha<- vec(t(u_alpha_2))
  
  all_u_delta<- c()
  
  for (j  in 1:(q-1) ){
    for (k in (j+1):q) {
      all_u_delta<- append(all_u_delta, u_delta_2[j,k,])
      
    }
    
  }
  all_u_theta<- c(all_u_alpha, all_u_delta)
  
  
  return(all_u_theta)
  
}


#MM algorithm for CMPLE
CMPLE.est<- function(mydata,epsilon,p,q){
  
  Covariates<- as.matrix(mydata[,1:p])
  Responses<- mydata[,(p+1):(p+q)]
  
  X<- model.matrix(~ Covariates)
  Y<- as.matrix(Responses)
  ##########
  #INITIALIZATION of model parameters#
  ################
  
  
  alpha<- array(0,dim=c(q,p+1))
  delta<- array(0,c(q,q,p+1))
  for (j in 1:q) {
    alpha[j,]<- rnorm(p+1,0,0.20)
    
  }
  
  for (j  in 1:(q-1) ){
    for (k in (j+1):q) {
      delta[j,k,]<- rnorm(p+1,0,0.20)
      
    }
    
  }
  
  all_alpha<- vec(t(alpha))
  
  all_delta<- c()
  
  for (j  in 1:(q-1) ){
    for (k in (j+1):q) {
      all_delta<- append(all_delta, delta[j,k,])
      
    }
    
  }
  theta<- c(all_alpha, all_delta)
  theta_est<-c()
  ###########
  print("Starting optimization. Please wait...")
  
  #####################
  diff<- 5
  
  count<-0
  while (abs(diff)>epsilon) {
    count<- count+1
    if(count>10000 ){
      break
    }
    #print(count)
    #print(diff)
    
    Sig<- array(0,dim=c(n,q))
    for (j in 1:q) {
      Sig[,j]<- exp (X %*% alpha[j,])
      
    }
    
    #Rho calculation 
    Rho<- array(0,c(n,q,q))
    
    for (j  in 1:(q-1) ){
      for (k in (j+1):q) {
        Rho[,j,k]<- 1- 2/(1+ exp(X %*% delta[j,k,]))
        
      }
      
    }
    
    for (i  in 1:n) {
      Rho[i,,]<- Rho[i,,]+ t(Rho[i,,])
      diag(Rho[i,,])<- 1
      
    }
    
    
    
    ######
    #first derivative calculation
    ######
    
    del_alpha<- array(0,c(n,q,q,p+1))
    
    
    for (i in 1:n) {
      
      for (j in 1: q)
      {
        for (k in 1:q) {
          if(k==j) next
          
          
          del_alpha[i,j,k,]<- (-1/Sig[i,j]+ (Y[i,j]**2/((Sig[i,j]**3)*(1-Rho[i,j,k]**2)))- (Y[i,j]*Y[i,k]/((Sig[i,j]**2)*(Sig[i,k])*(1-Rho[i,j,k]**2)))+
                                 (Y[i,j]*Y[i,k]/((Sig[i,j]**2)*(Sig[i,k])*(1+Rho[i,j,k]))))*Sig[i,j]* X[i,]
        }
        
      }
      
    }
    
    del_alpha2<- apply(del_alpha,c(2,4),sum)
    
    
    del_delta<- array(0,c(n,q,q,p+1))
    
    for(i in 1: n){
      for (j in 1: (q-1))
      {
        for (k in (j+1):q) {
          
          
          del_delta[i,j,k,]<-0.5*((Rho[i,j,k]/(1-Rho[i,j,k]**2))- ((Y[i,j]**2) * Rho[i,j,k]/(Sig[i,j]**2*(1-Rho[i,j,k]**2)**2)) -((Y[i,k]**2) * Rho[i,j,k]/((Sig[i,k]**2)*(1-Rho[i,j,k]**2)**2)) + (2*Y[i,j]*Y[i,k] * Rho[i,j,k]/(Sig[i,j]*Sig[i,k]*(1-Rho[i,j,k]**2)**2)) + (Y[i,j]*Y[i,k]/(Sig[i,j]*Sig[i,k]*(1+Rho[i,j,k])**2)))*(1-Rho[i,j,k]**2)*X[i,]
          
        }
      }
    }
    
    del_delta2<- apply(del_delta,c(2,3,4),sum)
    
    all_del_alpha<- vec(t(del_alpha2))
    all_del_delta<- c()
    
    for (j  in 1:(q-1) ){
      for (k in (j+1):q) {
        all_del_delta<- append(all_del_delta, del_delta2[j,k,])
        
      }
      
    }
    
    del_theta<- c(all_del_alpha,all_del_delta)
    
    ######
    #second derivative calculation
    ######
    del_del_alpha<- array(0,c(n,q,q,p+1,p+1))
    for (i in 1:n) {
      
      
      for (j in 1: q)
      {
        for (k in 1:q) {
          if(k==j) next
          del_del_alpha[i,j,k,,]<- ((-4* Y[i,j]**2/(Sig[i,j]**2*(1-Rho[i,j,k])**2)) - (3*(Y[i,j]**2+ Y[i,k]**2)/(2* Sig[i,j]* Sig[i,k]*(1-Rho[i,j,k])**2)) - (3*((Y[i,j] + Y[i,k])**2)/(2*Sig[i,j]*Sig[i,k]*(1+ Rho[i,j,k])))) *X[i,]%*%t(X[i,])
          
        }
      }
    }
    del_del_alpha2<- apply(del_del_alpha, c(2,4,5), sum)
    
    
    
    
    del_del_delta<- array(0,c(n,q,q,p+1,p+1))
    for (i in 1:n){
      
      for (j in 1: (q-1))
      {
        for (k in (j+1):q) {
          
          del_del_delta[i,j,k,,]<- (0.25*(((1+ Rho[i,j,k]**2)/(1-Rho[i,j,k]**2)**2) -(((Y[i,j]**2/Sig[i,j]**2)+ (Y[i,k]**2/Sig[i,k]**2))*((1+ 5* Rho[i,j,k]**2)/(1-Rho[i,j,k]**2)**3))+ ((1/(Sig[i,j]*Sig[i,k]*(1-Rho[i,j,k]**2)**3))*((Y[i,j]+Y[i,k])**2*(1+Rho[i,j,k]**2) -(Y[i,j]**2+ Y[i,k]**2)*(1+ 7* Rho[i,j,k]**2))) + ((Y[i,j]**2 + Y[i,k]**2- 4* (Y[i,j]+ Y[i,k])**2)/(2*Sig[i,j]*Sig[i,k]*(1+Rho[i,j,k])**3)))*((1-Rho[i,j,k]**2)**2)*X[i,]%*% t(X[i,]))
          - ((Rho[i,j,k]/(1-Rho[i,j,k]**2)- ((Y[i,j]**2) * Rho[i,j,k]/((Sig[i,j]**2)*(1-Rho[i,j,k]**2)**2)) -((Y[i,k]**2) * Rho[i,j,k]/((Sig[i,k]**2)*(1-Rho[i,j,k]**2)**2)) + 2*Y[i,j]*Y[i,k] * Rho[i,j,k]/(Sig[i,j]*Sig[i,k]*(1-Rho[i,j,k]**2)**2) + (Y[i,j]*Y[i,k])/(Sig[i,j]*Sig[i,k]*(1+Rho[i,j,k])**2)))*0.5*Rho[i,j,k]*(1-Rho[i,j,k]**2)*X[i,]%*%t(X[i,])
        }
      }
    }
    
    del_del_delta2<- apply(del_del_delta,c(2,3,4,5),sum)
    
    all_del_del_alpha<- c()
    for (i in 1:q) {
      all_del_del_alpha<- append(all_del_del_alpha,bdiag(del_del_alpha2[i,,]))
      
    }
    all_del_del_alpha<- as.matrix(bdiag(all_del_del_alpha))
    
    all_del_del_delta<- c()
    
    for (j  in 1:(q-1) ){
      for (k in (j+1):q) {
        all_del_del_delta<- append(all_del_del_delta,bdiag(del_del_delta2[j,k,,]))
        
      }
      
    }
    all_del_del_delta<- as.matrix(bdiag(all_del_del_delta))
    
    del_del_theta<- bdiag(all_del_del_alpha,all_del_del_delta)
    
    #####
    ##Update theta parameters
    #####
    increment<- solve(del_del_theta)%*% del_theta
    
    
    theta_new<- theta- increment
    
    old_theta<- theta
    old_theta[abs(theta)<0.03]<- 0.03
    
    
    diff<- sum(abs((theta-theta_new)/old_theta))   
    
    alpha_theta_new<- theta_new[1:(q*(p+1))]
    delta_theta_new<- theta_new[((q*(p+1))+1): length(theta_new)]
    
    place1=0
    for (i in 1:q) {
      
      alpha[i,]<- theta_new[(place1+1):(place1+p+1)]
      place1=place1+p+1
      
      
    }
    
    place2<- q*(p+1)
    for (j  in 1:(q-1) ){
      for (k in (j+1):q) {
        
        delta[j,k,]<- theta_new[(place2+1):(place2+p+1)]
        place2=place2+p+1
      }
    }
    
    theta<- theta_new
    
  }
  
  #####
  #store final estimates
  #####
  
  for (l in 1:(length(theta_new))) {
    theta_est[l]<- theta_new[l]
  }
  print("Parameter Estimates Done")
  return(theta_est)
  
}



CMPLE.std<- function(parameter,mydata,p,q){
  
  Covariates<- as.matrix(mydata[,1:p])
  Responses<- mydata[,(p+1):(p+q)]
  
  X<- model.matrix(~ Covariates)
  Y<- as.matrix(Responses)  
  
  n<- nrow(mydata)
  print("Calculating Standard Error of estimates. Please Wait...")
  u_new_theta<- jacobian(log_likelihood,parameter)
  print("Calculating jacobian. Please Wait...")
  
  J_theta<- (t(u_new_theta)%*% u_new_theta)/n 
  
  jaco<- (-jacobian(u_theta,parameter)/n)
  
  H_theta<- jaco
  
  G_inv_theta<- (solve(H_theta) %*% J_theta %*% solve(H_theta))/n 
  
  var_covar<- G_inv_theta
  std_error<- c()
  for (parms in 1:((q+1)*q*(p+1)/2)) {
    std_error[parms]<- (diag(var_covar)[parms])**0.5
    
  }
  print("Standard error calculation done")
  return(std_error)
  
}

CMPLE<- function(mydata,p,q,epsilon){
  suffix_alpha_1<- rep(1:q, each=p+1)
  suffix_alpha_2<- rep(0:p,q)
 # alpha_para<-paste0(expression(alpha),"_",suffix_alpha_1,",",suffix_alpha_2)
  alpha_para<-paste0(expression(alpha),"_",suffix_alpha_1,",",suffix_alpha_2)
  
  suffix_delta_1<-c()
  for (j  in 1:(q-1) ){
    for (k in (j+1):q) {
      suffix_delta_1<- append(suffix_delta_1,paste(j,k,sep=","))
    }
  }
  suffix_delta_1<- rep(suffix_delta_1,each=p+1)
  suffix_delta_2<- rep(0:p,q)
  delta_para<-paste0(expression(delta),"_",suffix_delta_1,",",suffix_delta_2)
  
  para.names<- c(alpha_para,delta_para)
  est<- CMPLE.est(mydata,epsilon,p,q)
  Estimates<- round(est,2)
  out1<- data.frame(para.names,Estimates)
  colnames(out1)<- c("Parameter","Estimates")
  std<- CMPLE.std(est, mydata,p,q)
  Standard_error<- round(std,2)
  out2<- data.frame(para.names,Standard_error)
  colnames(out2)<- c("Parameter","Standard Error")
  P.value<- round(2*(1-pnorm(abs(est/std))), 2)
  Lower_conf<- round(Estimates -1.96*Standard_error, 2)
  Upper_conf<- round(Estimates +1.96*Standard_error, 2)
  out3<- data.frame(para.names, Estimates, Standard_error, P.value, Lower_conf, Upper_conf)
  colnames(out3)<- c("Parameter","Estimates", "Standard_error", "P.value", "Lower_conf", "Upper_conf" )
  myout<- list(out1, out2, out3)
  names(myout)<- c("out1", "out2", "out3")
   return(myout)
  
}








