
#simulate data for 4 response and 10 predictors : 4 variance parameters(sigma):alpha : 4*11
#                               6 correlation parameters    (rho)  : 6* 11
#Total number of parameters     0.5*(p+1)*q*(q+1)
n= 1000
q<- 4 #number of phenotypes
p=10   #number of predictors 
set.seed(100)
x1<- rbinom(n,1,0.5)
x2<- rbinom(n,1,0.5)
x3<- rbinom(n,1,0.5)
x4<- rbinom(n,1,0.5)
x5<- rbinom(n,1,0.5)
x6<- rbinom(n,1,0.5)
x7<- rbinom(n,1,0.5)
x8<- rbinom(n,1,0.5)
x9<- rbinom(n,1,0.5)
x10<- rbinom(n,1,0.5)

#10 predictors 

x_mat<- model.matrix(~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)

#parameter values
#define alpha1, alpha2, alpha3, alpha4 (standard deviation parameters)


alpha1<- c(-1,1,0.20,-0.38,0.30,0,-0.50,-0.35,0.10,-0.40,1)
alpha2<- c(2, 0.20, -0.50, -0.41, 0, 0.35,0 ,0.30,0.19,-0.45,0)
alpha3<- c(1,0.30, 0.10, 0, -0.22, 0.15, -0.25, 1,-1,-0.15,0.25)
alpha4<- c(-1, -0.50, 0, 0.33,-0.18,-0.24, 0.18, 0, 0 ,0.55, 0.50)

#define delta12, delta13, delta14, delta23, delta24, delta34 (correlation parameters)

delta12<- c(0.25, 0.50, 0, -0.73, 0.30, 0, -0.83, 0.10,0, 0.15, 0.25)
delta13<- c(0.5, 0.25, 0.50, 0, 0.25,-0.55,0, -0.30, 0.20,0.45,-0.20)
delta14<- c(1, 0.25,-0.50,-0.25,0.50, 0.15,-0.15, -0.10, 0,-0.15,-0.25)
delta23<- c(0, 0, -0.25, 0, 0, 0, 0,-0.30, 0.20, 0, 0 )
delta24<- c(0.35, 0, -1, 0.45, -0.30, 0.50, 0.33, 0.10, -0.60, 0.55, -0.15)
delta34<- c(-1, 0.25,-0.25, -0.15, 0.23, 0.61,-0.15, 0,0.20, 0.30, -0.25)

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

mydata<- data.frame(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,response_df)