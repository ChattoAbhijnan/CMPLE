
#simulate data for 4 response : 4 variance parameters(sigma):alpha : 4*3
#                               6 correlation parameters    (rho)  : 6* 3
#Total number of parameters     0.5*(p+1)*q*(q+1)
n= 600
q<- 4 #number of phenotypes
p=2   #number of predictors 
set.seed(100)
x1<- rbinom(n,1,0.5)
x2<- rbinom(n,1,0.5)

#we have used 2 predictors for illustration purpose. 

x_mat<- model.matrix(~x1+x2)

#parameter values
#define alpha1, alpha2, alpha3, alpha4 (standard deviation parameters)


alpha1<- c(-1,1,0.20)
alpha2<- c(2,0.20,-0.50)
alpha3<- c(1,0.30,0.10)
alpha4<- c(-1,-0.50,1)

#define delta12, delta13, delta14, delta23, delta24, delta34 (correlation parameters)

delta12<- c(0.25,0.50,1)
delta13<- c(0.25,0.25,0.50)
delta14<- c(0.25,0.25,-0.50)
delta23<- c(-0.05,1,-0.25)
delta24<- c(0.25,0.50,-1)
delta34<- c(-0.05,0.25,-0.25)

sigma1<- exp(x_mat %*% alpha1)
sigma2<- exp(x_mat %*% alpha2)
sigma3<- exp(x_mat %*% alpha3)
sigma4<- exp(x_mat %*% alpha4)

rho12<- 1- 2/(1+exp(x_mat %*% delta12))
rho13<- 1- 2/(1+exp(x_mat %*% delta13))
rho14<- 1- 2/(1+exp(x_mat %*% delta14))
rho23<- 1- 2/(1+exp(x_mat %*% delta23))
rho24<- 1- 2/(1+exp(x_mat %*% delta24))
rho34<- 1- 2/(1+exp(x_mat %*% delta34))
####

####

y<- list()
for (s in 1:n) {
  A<- array(0,dim=c(q,q))
  diag(A)<- c(sigma1[s],sigma2[s],sigma3[s],sigma4[s])
  R<- array(0,dim=c(q,q))
  R[upper.tri(R,diag = FALSE)]<-
    c(rho12[s],rho13[s],rho23[s],rho14[s],rho24[s],rho34[s])
  rho<- R + t(R)
  diag(rho)<- rep(1,q)
  rho_1<-nearPD(rho, corr = TRUE, keepDiag = TRUE)
  rho_1<- as.matrix(rho_1[[1]])
  Omega<- A %*% rho_1 %*% A
  mu<- c(0,0,0,0)
  y[[s]]<-mvrnorm(n = 1, mu,  Omega)
}
response_df<- data.frame(matrix(unlist(y), nrow=length(y), byrow=TRUE))
colnames(response_df)<- c("Y1","Y2","Y3","Y4")

Y<- as.matrix(response_df)
X<- x_mat

mydata<- data.frame(x1,x2,response_df)
